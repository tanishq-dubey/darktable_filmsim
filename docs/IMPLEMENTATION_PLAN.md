# darktable Film Simulation Plugin - Complete Implementation Plan

This document provides a step-by-step guide to implementing the full-featured darktable native plugin.

## Architecture Summary

**Goal**: Port cppfilmsim's complete film simulation pipeline to a darktable IOP module that outputs film negative density (for use with negadoctor).

**Pipeline Flow**:
```
Linear Scene RGB → Halation → Film Response (H-D) → DIR → Dye Mix → Negative Density
```

**Output**: Film negative in dye density space (D_C, D_M, D_Y) stored in RGB channels

## Implementation Phases

### Phase 1: Core Infrastructure ✓ STARTED

**Files Created**:
- `README.md` - User documentation
- `filmsim_common.h` - C data structures
- `IMPLEMENTATION_PLAN.md` - This file

**Remaining**:
- `CMakeLists.txt` - Build configuration
- `filmsim_curves.c` - Monotone cubic spline (port from `Spline.h`)
- `filmsim_data.c` - Film profile loading (JSON parsing)

### Phase 2: Effects Implementation

#### 2a. Monotone Cubic Spline (`filmsim_curves.c`)

**Source**: `src/Spline.h` (C++ template class)

**Port Strategy**:
```c
// filmsim_curves.c

#include "filmsim_common.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

void filmsim_spline_init(filmsim_spline_t *spline) {
    spline->x = NULL;
    spline->y = NULL;
    spline->m = NULL;
    spline->n = 0;
}

void filmsim_spline_fit(filmsim_spline_t *spline, 
                         const double *x, const double *y, int n) {
    // Allocate memory
    spline->n = n;
    spline->x = (double*)malloc(n * sizeof(double));
    spline->y = (double*)malloc(n * sizeof(double));
    spline->m = (double*)malloc(n * sizeof(double));
    
    // Copy data
    memcpy(spline->x, x, n * sizeof(double));
    memcpy(spline->y, y, n * sizeof(double));
    
    // Compute secants
    double *d = (double*)malloc((n-1) * sizeof(double));
    for (int i = 0; i < n-1; i++) {
        double h = spline->x[i+1] - spline->x[i];
        d[i] = (spline->y[i+1] - spline->y[i]) / h;
    }
    
    // Initialize tangents (Fritsch-Carlson algorithm)
    spline->m[0] = d[0];
    spline->m[n-1] = d[n-2];
    
    for (int i = 1; i < n-1; i++) {
        if (d[i-1] * d[i] <= 0) {
            spline->m[i] = 0;
        } else {
            double h1 = spline->x[i] - spline->x[i-1];
            double h2 = spline->x[i+1] - spline->x[i];
            double w1 = 2*h2 + h1;
            double w2 = h2 + 2*h1;
            spline->m[i] = (w1 + w2) / (w1/d[i-1] + w2/d[i]);
        }
    }
    
    // Ensure monotonicity
    for (int i = 0; i < n-1; i++) {
        if (d[i] == 0) {
            spline->m[i] = 0;
            spline->m[i+1] = 0;
        } else {
            double alpha = spline->m[i] / d[i];
            double beta = spline->m[i+1] / d[i];
            double s = alpha*alpha + beta*beta;
            if (s > 9) {
                double tau = 3 / sqrt(s);
                spline->m[i] = tau * alpha * d[i];
                spline->m[i+1] = tau * beta * d[i];
            }
        }
    }
    
    free(d);
}

double filmsim_spline_interpolate(const filmsim_spline_t *spline, double x) {
    if (spline->n == 0) return 0.0;
    if (x <= spline->x[0]) return spline->y[0];
    if (x >= spline->x[spline->n-1]) return spline->y[spline->n-1];
    
    // Binary search for segment
    int i = 0;
    for (int j = spline->n-1; i < j; ) {
        int m = (i + j) / 2;
        if (x < spline->x[m]) j = m;
        else i = m + 1;
    }
    i--;
    
    // Hermite interpolation
    double h = spline->x[i+1] - spline->x[i];
    double t = (x - spline->x[i]) / h;
    double t2 = t*t;
    double t3 = t2*t;
    
    double h00 = 2*t3 - 3*t2 + 1;
    double h10 = t3 - 2*t2 + t;
    double h01 = -2*t3 + 3*t2;
    double h11 = t3 - t2;
    
    return h00*spline->y[i] + h10*h*spline->m[i] + 
           h01*spline->y[i+1] + h11*h*spline->m[i+1];
}

void filmsim_spline_free(filmsim_spline_t *spline) {
    free(spline->x);
    free(spline->y);
    free(spline->m);
    spline->x = spline->y = spline->m = NULL;
    spline->n = 0;
}
```

#### 2b. Halation Effect (`filmsim_halation.c`)

**Source**: `src/Halation.cpp` (341 lines, already well-structured C++)

**Port Strategy**:
- Replace `std::vector<float>` with `float*` arrays
- Replace `ScopedTimer` with manual timing (optional in plugin)
- Keep FFTW code identical (C library, no changes needed)
- Replace `#pragma omp parallel for` → same (OpenMP is C-compatible)
- Remove `std::cout` → use darktable's `dt_print()` or silence

**Key Changes**:
```c
// Instead of:
void applyHalation(std::vector<float>& imageData, ...)

// Use:
void filmsim_apply_halation(float *imageData, int width, int height, ...)
```

**Critical**: Keep downsample factor (DS=4), FFT logic, bilinear upsampling exactly as-is.

#### 2c. DIR Effect (`filmsim_dir.c`)

**Source**: `src/DIR.cpp` (364 lines, similar structure to Halation)

**Port Strategy**: Same as halation - mostly mechanical translation from C++ to C.

**Key Points**:
- Downsample factor DS=2 (vs 4 for halation)
- Gaussian PSF (vs exponential for halation)
- Same FFTW convolution approach
- DIR formula: `adjusted = original + amount * (original - blurred)`

### Phase 3: Film Profile System

#### 3a. JSON Parsing (`filmsim_data.c`)

**Challenge**: cppfilmsim uses nlohmann/json (C++). darktable doesn't bundle this.

**Solution Options**:

**Option A: GLib JSON Parser** (Recommended - already in darktable)
```c
#include <glib.h>
#include <json-glib/json-glib.h>

bool filmsim_load_profile_from_json(filmsim_film_profile_t *profile, 
                                     const char *json_string) {
    JsonParser *parser = json_parser_new();
    GError *error = NULL;
    
    if (!json_parser_load_from_data(parser, json_string, -1, &error)) {
        g_error_free(error);
        g_object_unref(parser);
        return false;
    }
    
    JsonNode *root = json_parser_get_root(parser);
    JsonObject *obj = json_node_get_object(root);
    
    // Parse "properties" -> "info" -> "name"
    JsonObject *props = json_object_get_object_member(obj, "properties");
    JsonObject *info = json_object_get_object_member(props, "info");
    const char *name = json_object_get_string_member(info, "name");
    strncpy(profile->info.name, name, sizeof(profile->info.name)-1);
    
    // ... parse all other fields ...
    
    g_object_unref(parser);
    return true;
}
```

**Option B: Pre-compile JSON to C** (Build-time conversion)
```python
# Python script: json_to_c.py
import json
films = []
for json_file in glob("../data/*.json"):
    with open(json_file) as f:
        data = json.load(f)
        # Generate C initializer
        print(f"static const filmsim_film_profile_t film_{name} = {{")
        print(f"  .info.name = \"{data['properties']['info']['name']}\",")
        # ... etc
        print("};")
```

**Recommendation**: Use Option A initially (GLib JSON), consider Option B for performance later.

#### 3b. Embedded Film Profiles

**Current Setup**: cppfilmsim uses CMake Resource Compiler (cmrc) to embed `data/*.json`

**darktable Approach**: Same! darktable also uses cmrc for resources.

**CMakeLists.txt** (partial):
```cmake
# Embed all film profiles
file(GLOB FILM_PROFILES "${CMAKE_SOURCE_DIR}/data/*.json")

cmrc_add_resource_library(
    filmsim_profiles
    WHENCE ${CMAKE_SOURCE_DIR}/data
    ${FILM_PROFILES}
)

# Link to plugin
add_library(filmsim MODULE filmsim.c ...)
target_link_libraries(filmsim cmrc::filmsim_profiles)
```

**Access in C**:
```c
#include <cmrc/cmrc.h>

CMRC_DECLARE(filmsim_profiles);

const char* filmsim_get_embedded_profile(const char *film_name) {
    auto fs = cmrc::filmsim_profiles::get_filesystem();
    char path[256];
    snprintf(path, sizeof(path), "%s.json", film_name);
    
    auto file = fs.open(path);
    return file.begin();  // Returns const char* to embedded data
}
```

### Phase 4: Main IOP Module (`filmsim.c`)

This is the darktable-specific integration layer. It implements darktable's IOP API.

**Structure** (simplified, ~500-800 lines total):

```c
#include "filmsim_common.h"
#include "develop/imageop.h"
#include "bauhaus/bauhaus.h"
#include "gui/gtk.h"

/* Module metadata */
DT_MODULE(1)  // Version 1

const char *name() {
    return _("film simulation");
}

const char *aliases() {
    return _("filmsim|analog|negative");
}

const char *description() {
    return _("Simulate analog film response with halation, DIR, and dye mixing");
}

int default_group() {
    return IOP_GROUP_COLOR | IOP_GROUP_EFFECT;
}

int flags() {
    return IOP_FLAGS_SUPPORTS_BLENDING | IOP_FLAGS_ALLOW_TILING;
}

dt_iop_colorspace_type_t default_colorspace() {
    return IOP_CS_RGB;  // Linear RGB
}

/* Parameters (stored in database/XMP) */
typedef struct dt_iop_filmsim_params_t {
    int film_profile;           // Index into embedded films (0-17)
    float push_pull;            // -3.0 to +3.0
    int enable_halation;        // 0 or 1
    int enable_dir;             // 0 or 1 (default 1)
    int enable_dye_mix;         // 0 or 1 (default 1)
    float pixel_size;           // 1.0 to 20.0 (default 6.0)
    // ... more advanced params
} dt_iop_filmsim_params_t;

/* GUI data (runtime only) */
typedef struct dt_iop_filmsim_gui_data_t {
    GtkWidget *film_combo;
    GtkWidget *push_pull_slider;
    GtkWidget *halation_toggle;
    GtkWidget *dir_toggle;
    GtkWidget *dye_mix_toggle;
    GtkWidget *pixel_size_slider;
} dt_iop_filmsim_gui_data_t;

/* Processing data (cached, per-pipeline instance) */
typedef struct dt_iop_filmsim_data_t {
    filmsim_film_profile_t film;  // Loaded film profile
    float push_pull;
    int enable_halation;
    int enable_dir;
    int enable_dye_mix;
    float pixel_size;
    
    // Precomputed LUTs for performance
    float film_lut_r[65536];
    float film_lut_g[65536];
    float film_lut_b[65536];
} dt_iop_filmsim_data_t;

/* MAIN PROCESSING FUNCTION */
void process(struct dt_iop_module_t *self,
             dt_dev_pixelpipe_iop_t *piece,
             const void *const ivoid,
             void *const ovoid,
             const dt_iop_roi_t *const roi_in,
             const dt_iop_roi_t *const roi_out)
{
    dt_iop_filmsim_data_t *data = (dt_iop_filmsim_data_t*)piece->data;
    
    const float *in = (const float*)ivoid;
    float *out = (float*)ovoid;
    
    const int ch = piece->colors;  // 4 (RGBA)
    const int width = roi_out->width;
    const int height = roi_out->height;
    const size_t npixels = (size_t)width * height;
    
    // Copy input to output (working buffer)
    memcpy(out, in, sizeof(float) * ch * npixels);
    
    // 1. Apply halation (optional, in exposure domain)
    if (data->enable_halation) {
        filmsim_apply_halation(out, width, height, 
                               &data->film.halation, data->pixel_size);
    }
    
    // 2. Apply film response curves (H-D) → density
    //    Uses precomputed LUTs for speed
    DT_OMP_FOR()
    for (size_t k = 0; k < npixels; k++) {
        float r = out[k*ch + 0];
        float g = out[k*ch + 1];
        float b = out[k*ch + 2];
        
        // Clamp and quantize to LUT index [0, 65535]
        int idx_r = CLAMP((int)(r * 65535.0f), 0, 65535);
        int idx_g = CLAMP((int)(g * 65535.0f), 0, 65535);
        int idx_b = CLAMP((int)(b * 65535.0f), 0, 65535);
        
        // Lookup density
        out[k*ch + 0] = data->film_lut_r[idx_r];  // D_C (cyan)
        out[k*ch + 1] = data->film_lut_g[idx_g];  // D_M (magenta)
        out[k*ch + 2] = data->film_lut_b[idx_b];  // D_Y (yellow)
        // Alpha unchanged
    }
    
    // 3. Apply DIR (optional, in density domain)
    if (data->enable_dir) {
        filmsim_apply_dir(out, width, height, 
                          &data->film.couplers, data->pixel_size);
    }
    
    // 4. Apply dye mixing (optional)
    if (data->enable_dye_mix) {
        filmsim_apply_dye_mixing(out, width, height, &data->film.couplers);
    }
    
    // Output is now film negative density (D_C, D_M, D_Y)
    // User applies negadoctor module next to invert to positive
}

/* Precompute film curve LUTs */
void commit_params(struct dt_iop_module_t *self, ...) {
    dt_iop_filmsim_data_t *data = piece->data;
    dt_iop_filmsim_params_t *p = params;
    
    // Load film profile
    const char *film_json = filmsim_get_embedded_profile(p->film_profile);
    filmsim_load_profile_from_json(&data->film, film_json);
    
    // Copy params
    data->push_pull = p->push_pull;
    data->enable_halation = p->enable_halation;
    // ... etc
    
    // Build LUT from splines
    filmsim_spline_t spline_r, spline_g, spline_b;
    // ... fit splines from film curve data ...
    
    float push_pull_log = p->push_pull * 0.301f;  // 1 stop = 0.301 log units
    
    for (int i = 0; i < 65536; i++) {
        float exposure = (float)i / 65535.0f;
        float logE = log10f(exposure + 1e-8f) + push_pull_log;
        
        data->film_lut_r[i] = filmsim_spline_interpolate(&spline_r, logE);
        data->film_lut_g[i] = filmsim_spline_interpolate(&spline_g, logE);
        data->film_lut_b[i] = filmsim_spline_interpolate(&spline_b, logE);
    }
    
    filmsim_spline_free(&spline_r);
    // ... etc
}

/* GUI initialization */
void gui_init(struct dt_iop_module_t *self) {
    dt_iop_filmsim_gui_data_t *g = IOP_GUI_ALLOC(filmsim);
    
    // Film selection combo box
    g->film_combo = dt_bauhaus_combobox_new(self);
    dt_bauhaus_widget_set_label(g->film_combo, NULL, _("film profile"));
    
    // Add all 18 films
    dt_bauhaus_combobox_add(g->film_combo, "Kodak Portra 400");
    dt_bauhaus_combobox_add(g->film_combo, "Kodak Ektar 100");
    dt_bauhaus_combobox_add(g->film_combo, "Fuji Pro H 400");
    // ... all 18 films ...
    
    g_signal_connect(G_OBJECT(g->film_combo), "value-changed",
                     G_CALLBACK(film_changed), self);
    
    gtk_box_pack_start(GTK_BOX(self->widget), g->film_combo, TRUE, TRUE, 0);
    
    // Push/pull slider
    g->push_pull_slider = dt_bauhaus_slider_new_with_range(self, -3.0, 3.0, 0.1, 0.0, 2);
    dt_bauhaus_widget_set_label(g->push_pull_slider, NULL, _("push/pull"));
    dt_bauhaus_slider_set_format(g->push_pull_slider, _("%.1f stops"));
    g_signal_connect(G_OBJECT(g->push_pull_slider), "value-changed",
                     G_CALLBACK(push_pull_changed), self);
    gtk_box_pack_start(GTK_BOX(self->widget), g->push_pull_slider, TRUE, TRUE, 0);
    
    // Halation toggle
    g->halation_toggle = dt_bauhaus_toggle_new(self);
    dt_bauhaus_widget_set_label(g->halation_toggle, NULL, _("halation"));
    g_signal_connect(G_OBJECT(g->halation_toggle), "value-changed",
                     G_CALLBACK(halation_changed), self);
    gtk_box_pack_start(GTK_BOX(self->widget), g->halation_toggle, TRUE, TRUE, 0);
    
    // ... DIR toggle, dye mix toggle, pixel size slider ...
}

/* Callbacks for GUI changes */
static void film_changed(GtkWidget *widget, dt_iop_module_t *self) {
    dt_iop_filmsim_params_t *p = self->params;
    p->film_profile = dt_bauhaus_combobox_get(widget);
    dt_dev_add_history_item(darktable.develop, self, TRUE);
}

// ... similar for other controls ...
```

### Phase 5: Build System

**CMakeLists.txt** (external plugin version):

```cmake
cmake_minimum_required(VERSION 3.10)
project(filmsim_plugin VERSION 1.0)

set(CMAKE_C_STANDARD 11)

# Find dependencies
find_package(PkgConfig REQUIRED)
pkg_check_modules(GLIB REQUIRED glib-2.0)
pkg_check_modules(GTK REQUIRED gtk+-3.0)
pkg_check_modules(JSON_GLIB REQUIRED json-glib-1.0)

find_package(FFTW3f REQUIRED)
find_package(OpenMP REQUIRED)

# darktable SDK (user must set DARKTABLE_INCLUDE_DIR)
if(NOT DEFINED DARKTABLE_INCLUDE_DIR)
    message(FATAL_ERROR "Please set -DDARKTABLE_INCLUDE_DIR=/path/to/darktable/src")
endif()

# Source files
set(PLUGIN_SOURCES
    filmsim.c
    filmsim_curves.c
    filmsim_halation.c
    filmsim_dir.c
    filmsim_data.c
)

# Embed film profiles
include(FetchContent)
FetchContent_Declare(
    cmrc
    GIT_REPOSITORY https://github.com/vector-of-bool/cmrc.git
    GIT_TAG 2.0.1
)
FetchContent_MakeAvailable(cmrc)

file(GLOB FILM_PROFILES "${CMAKE_SOURCE_DIR}/../data/*.json")
cmrc_add_resource_library(
    filmsim_profiles
    WHENCE ${CMAKE_SOURCE_DIR}/../data
    ${FILM_PROFILES}
)

# Build plugin as shared library
add_library(filmsim MODULE ${PLUGIN_SOURCES})

target_include_directories(filmsim PRIVATE
    ${DARKTABLE_INCLUDE_DIR}
    ${GLIB_INCLUDE_DIRS}
    ${GTK_INCLUDE_DIRS}
    ${JSON_GLIB_INCLUDE_DIRS}
)

target_link_libraries(filmsim PRIVATE
    ${GLIB_LIBRARIES}
    ${GTK_LIBRARIES}
    ${JSON_GLIB_LIBRARIES}
    FFTW3::fftw3f
    OpenMP::OpenMP_C
    cmrc::filmsim_profiles
)

# Install to darktable plugins directory
install(TARGETS filmsim 
    LIBRARY DESTINATION $ENV{HOME}/.config/darktable/plugins
)
```

**Build & Install**:
```bash
cd darktable_plugin
mkdir build && cd build
cmake .. -DDARKTABLE_INCLUDE_DIR=/path/to/darktable/src
make -j
make install  # Copies to ~/.config/darktable/plugins/
```

## Testing Workflow

1. **Build plugin**: `make && make install`
2. **Launch darktable**: `darktable -d opencl,perf` (for performance monitoring)
3. **Load RAW file**
4. **Enable modules**:
   - Disable auto-applied modules if needed
   - Enable: `input color profile → filmsim → negadoctor`
5. **Configure**:
   - filmsim: Select film (e.g., "Kodak Portra 400")
   - filmsim: Enable halation, DIR
   - negadoctor: Set D-min, scanner, etc.
6. **Compare** output to standalone cppfilmsim

## Performance Expectations

- **Film curve LUT**: ~100 µs (one-time in commit_params)
- **Halation** (4K image, DS=4): ~50-150 ms (FFTW-dependent)
- **DIR** (4K image, DS=2): ~100-300 ms
- **Dye mixing**: ~10 ms (simple matrix multiply)
- **Total pipeline**: 200-500 ms for 4K image (acceptable for interactive editing)

## Known Challenges

1. **FFTW Thread Safety**: Ensure `fftwf_init_threads()` called in `init_global()`
2. **Memory Management**: No `std::vector` - must manually malloc/free
3. **Tiling Support**: May need `modify_roi_in()` for halation/DIR kernels
4. **Color Space**: darktable's RGB may differ from cppfilmsim's assumed sRGB primaries

## Next Steps After Implementation

1. **Upstream to darktable**: Submit PR to darktable/darktable repo
2. **Documentation**: Write user manual section
3. **Presets**: Create presets for common film + development combinations
4. **GPU/OpenCL**: Port halation/DIR to OpenCL for GPU acceleration

## File Checklist

- [x] `README.md`
- [x] `IMPLEMENTATION_PLAN.md`
- [x] `filmsim_common.h`
- [ ] `CMakeLists.txt`
- [ ] `filmsim_curves.c`
- [ ] `filmsim_halation.c`
- [ ] `filmsim_dir.c`
- [ ] `filmsim_data.c`
- [ ] `filmsim.c` (main IOP module)

**Estimated Implementation Time**: 20-30 hours for experienced C/darktable developer

## Questions/Decisions Needed

None remaining - ready to implement!
