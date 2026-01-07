/*
    This file is part of darktable,
    Film Simulation IOP module

    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "bauhaus/bauhaus.h"
#include "common/darktable.h"
#include "common/imagebuf.h"
#include "control/control.h"
#include "develop/develop.h"
#include "develop/imageop.h"
#include "develop/imageop_math.h"
#include "develop/imageop_gui.h"
#include "develop/pixelpipe.h"
#include "gui/gtk.h"
#include "iop/iop_api.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>

#include "filmsim_common.h"
#include "filmsim_spectral.h"

DT_MODULE_INTROSPECTION(3, dt_iop_filmsim_params_t)

#define FILMSIM_LUT_SIZE 65536
#define FILMSIM_NUM_FILMS 18

static const char *filmsim_film_names[] = {
    "kodak_portra_400",
    "kodak_ektar_100",
    "kodak_ultramax_400",
    "kodak_vision3_500t",
    "kodak_kodachrome_64",
    "kodak_ektachrome_p_1600",
    "kodak_royal_gold_1000",
    "kodak_portra_100T",
    "kodak_t-max_400_tmax-rs_7min",
    "fuji_proh_400",
    "fuji_superia_1600",
    "fuji_superia_reala_100",
    "fuji_fujichrome_velvia_100",
    "fuji_fujichrome_astia_f_100",
    "fuji_neopan_arcos_ii_100",
    "cinestill_800t",
    "agfacolor_optima_ii_100",
    "agfacolor_portrait_xps_160",
    NULL
};

typedef struct dt_iop_filmsim_params_t
{
    int32_t film_profile;       // $DEFAULT: 0 $DESCRIPTION: "film profile"
    float push_pull;            // $MIN: -10.0 $MAX: 10.0 $DEFAULT: 0.0 $DESCRIPTION: "push/pull"
    gboolean enable_halation;   // $DEFAULT: TRUE $DESCRIPTION: "halation"
    gboolean enable_dir;        // $DEFAULT: TRUE $DESCRIPTION: "DIR"
    gboolean enable_dye_mix;    // $DEFAULT: TRUE $DESCRIPTION: "dye mixing"
    gboolean enable_spectral;   // $DEFAULT: FALSE $DESCRIPTION: "spectral mode"
    float halation_strength;    // $MIN: 0.0 $MAX: 2.0 $DEFAULT: 1.0 $DESCRIPTION: "halation strength"
    float dir_strength;         // $MIN: 0.0 $MAX: 2.0 $DEFAULT: 1.0 $DESCRIPTION: "DIR strength"
    float pixel_size;           // $MIN: 1.0 $MAX: 20.0 $DEFAULT: 6.0 $DESCRIPTION: "pixel size"
} dt_iop_filmsim_params_t;

typedef struct dt_iop_filmsim_gui_data_t
{
    GtkWidget *film_combo;
    GtkWidget *push_pull_slider;
    GtkWidget *halation_toggle;
    GtkWidget *dir_toggle;
    GtkWidget *dye_mix_toggle;
    GtkWidget *spectral_toggle;
    GtkWidget *halation_strength_slider;
    GtkWidget *dir_strength_slider;
    GtkWidget *pixel_size_slider;
} dt_iop_filmsim_gui_data_t;

typedef struct dt_iop_filmsim_data_t
{
    filmsim_film_profile_t film;
    float push_pull;
    gboolean enable_halation;
    gboolean enable_dir;
    gboolean enable_dye_mix;
    gboolean enable_spectral;
    gboolean has_spectral_data;  /* Whether film profile has spectral data */
    float halation_strength;
    float dir_strength;
    float pixel_size;
    float film_lut_r[FILMSIM_LUT_SIZE];
    float film_lut_g[FILMSIM_LUT_SIZE];
    float film_lut_b[FILMSIM_LUT_SIZE];
    filmsim_spectral_lut_t spectral_lut;
} dt_iop_filmsim_data_t;

/* Global RGB-to-spectral model (loaded once, shared by all instances) */
static filmsim_rgb2spec_t g_rgb2spec_model = {0};
static gboolean g_rgb2spec_loaded = FALSE;

const char *name()
{
    return _("film simulation");
}

const char *aliases()
{
    return _("filmsim|analog|negative|halation");
}

const char** description(dt_iop_module_t *self)
{
    return dt_iop_set_description(self,
        _("simulate analog film with H-D curves, halation, and DIR effects"),
        _("creative"),
        _("linear, RGB, scene-referred"),
        _("linear, RGB"),
        _("linear, RGB, scene-referred"));
}

int default_group()
{
    return IOP_GROUP_EFFECT;
}

int flags()
{
    return IOP_FLAGS_SUPPORTS_BLENDING;
}

dt_iop_colorspace_type_t default_colorspace(dt_iop_module_t *self,
                                             dt_dev_pixelpipe_t *pipe,
                                             dt_dev_pixelpipe_iop_t *piece)
{
    return IOP_CS_RGB;
}

static void build_film_lut(dt_iop_filmsim_data_t *data)
{
    filmsim_film_profile_t *film = &data->film;
    
    if(film->num_curves < 2)
    {
        /* No valid curves - use identity mapping */
        for(int i = 0; i < FILMSIM_LUT_SIZE; i++)
        {
            float v = (float)i / (FILMSIM_LUT_SIZE - 1);
            data->film_lut_r[i] = v;
            data->film_lut_g[i] = v;
            data->film_lut_b[i] = v;
        }
        return;
    }

    int n = film->num_curves;
    double *x_vals = (double *)malloc(n * sizeof(double));
    double *y_r = (double *)malloc(n * sizeof(double));
    double *y_g = (double *)malloc(n * sizeof(double));
    double *y_b = (double *)malloc(n * sizeof(double));

    for(int i = 0; i < n; i++)
    {
        x_vals[i] = film->curves.rgb[i].d;
        y_r[i] = film->curves.rgb[i].r;
        y_g[i] = film->curves.rgb[i].g;
        y_b[i] = film->curves.rgb[i].b;
    }

    filmsim_spline_t spline_r, spline_g, spline_b;
    filmsim_spline_fit(&spline_r, x_vals, y_r, n);
    filmsim_spline_fit(&spline_g, x_vals, y_g, n);
    filmsim_spline_fit(&spline_b, x_vals, y_b, n);

    /* Calculate exposure scale k from calibration
     * middle_gray_logE corresponds to 0.18 linear value
     * log10(0.18 * k) = middle_gray_logE
     * k = 10^(middle_gray_logE) / 0.18
     */
    const double LN10 = 2.302585092994046;
    double k = exp(film->calibration.middle_gray_logE * LN10) / 0.18;
    
    /* Push/Pull shift: 1 stop = log10(2) = 0.301 */
    double delta_H = data->push_pull * log10(2.0);

    for(int i = 0; i < FILMSIM_LUT_SIZE; i++)
    {
        double linear = (double)i / (FILMSIM_LUT_SIZE - 1);
        double val = fmax(linear, 1e-6);
        
        /* Convert linear RGB to log exposure using calibration scale */
        double E = val * k;
        double logE = log10(E) + delta_H;

        data->film_lut_r[i] = (float)filmsim_spline_interpolate(&spline_r, logE);
        data->film_lut_g[i] = (float)filmsim_spline_interpolate(&spline_g, logE);
        data->film_lut_b[i] = (float)filmsim_spline_interpolate(&spline_b, logE);
    }

    filmsim_spline_free(&spline_r);
    filmsim_spline_free(&spline_g);
    filmsim_spline_free(&spline_b);
    free(x_vals);
    free(y_r);
    free(y_g);
    free(y_b);
}

void process(struct dt_iop_module_t *self,
             dt_dev_pixelpipe_iop_t *piece,
             const void *const ivoid,
             void *const ovoid,
             const dt_iop_roi_t *const roi_in,
             const dt_iop_roi_t *const roi_out)
{
    dt_iop_filmsim_data_t *data = (dt_iop_filmsim_data_t *)piece->data;

    const int ch = piece->colors;
    const int width = roi_out->width;
    const int height = roi_out->height;
    const size_t npixels = (size_t)width * height;

    const float *in = (const float *)ivoid;
    float *out = (float *)ovoid;

    /* Create 3-channel working buffer for effects that need it */
    float *work = dt_alloc_align_float(npixels * 3);
    if(!work)
    {
        /* Fallback: just copy input to output */
        dt_iop_copy_image_roi(ovoid, ivoid, ch, roi_in, roi_out);
        return;
    }

    /* Extract RGB from RGBA into working buffer */
    DT_OMP_FOR()
    for(size_t k = 0; k < npixels; k++)
    {
        work[k * 3 + 0] = in[k * ch + 0];
        work[k * 3 + 1] = in[k * ch + 1];
        work[k * 3 + 2] = in[k * ch + 2];
    }

    /* Step 0 (optional): Apply spectral exposure conversion */
    if(data->enable_spectral && data->spectral_lut.valid)
    {
        filmsim_spectral_apply_exposure(&data->spectral_lut, work, width, height);
    }

    /* Step 1: Apply halation (in exposure/linear domain, before H-D curves) */
    if(data->enable_halation && data->halation_strength > 0.001f)
    {
        filmsim_halation_t scaled_halation = data->film.halation;
        scaled_halation.strength.r *= data->halation_strength;
        scaled_halation.strength.g *= data->halation_strength;
        scaled_halation.strength.b *= data->halation_strength;
        
        filmsim_apply_halation(work, width, height, &scaled_halation, data->pixel_size);
    }

    /* Step 2: Apply H-D curves (exposure -> density) */
    const gboolean is_mono = (data->film.info.film_type == FILMSIM_FILM_TYPE_MONOCHROME_NEGATIVE);
    
    DT_OMP_FOR()
    for(size_t k = 0; k < npixels; k++)
    {
        float r = work[k * 3 + 0];
        float g = work[k * 3 + 1];
        float b = work[k * 3 + 2];

        if(is_mono)
        {
            /* B&W film: convert RGB to mono luminance first, then apply curve */
            /* Using Rec. 709 weights - could be improved with spectral sensitivity */
            float mono = r * 0.2126f + g * 0.7152f + b * 0.0722f;
            int idx = CLAMP((int)(mono * (FILMSIM_LUT_SIZE - 1)), 0, FILMSIM_LUT_SIZE - 1);
            float result = data->film_lut_r[idx];  /* All LUTs identical for mono */
            
            work[k * 3 + 0] = result;
            work[k * 3 + 1] = result;
            work[k * 3 + 2] = result;
        }
        else
        {
            /* Color film: apply curves per channel */
            int idx_r = CLAMP((int)(r * (FILMSIM_LUT_SIZE - 1)), 0, FILMSIM_LUT_SIZE - 1);
            int idx_g = CLAMP((int)(g * (FILMSIM_LUT_SIZE - 1)), 0, FILMSIM_LUT_SIZE - 1);
            int idx_b = CLAMP((int)(b * (FILMSIM_LUT_SIZE - 1)), 0, FILMSIM_LUT_SIZE - 1);

            work[k * 3 + 0] = data->film_lut_r[idx_r];
            work[k * 3 + 1] = data->film_lut_g[idx_g];
            work[k * 3 + 2] = data->film_lut_b[idx_b];
        }
    }

    /* Step 3: Apply DIR (in density domain, after H-D curves) */
    if(data->enable_dir && data->dir_strength > 0.001f)
    {
        filmsim_couplers_t scaled_couplers = data->film.couplers;
        scaled_couplers.dir_amount_rgb[0] *= data->dir_strength;
        scaled_couplers.dir_amount_rgb[1] *= data->dir_strength;
        scaled_couplers.dir_amount_rgb[2] *= data->dir_strength;
        
        filmsim_apply_dir(work, width, height, &scaled_couplers, data->pixel_size);
    }

    /* Step 4: Apply dye mixing (in density domain, color films only) */
    if(data->enable_dye_mix && !is_mono)
    {
        filmsim_apply_dye_mixing(work, width, height, &data->film.couplers);
    }

    /* Step 5: Convert density to transmittance (or use spectral dye model) */
    if(data->enable_spectral && data->spectral_lut.valid)
    {
        /* Use spectral dye model for more accurate output */
        filmsim_spectral_apply_dye(&data->spectral_lut, work, width, height);
        
        /* Copy to output with alpha */
        DT_OMP_FOR()
        for(size_t k = 0; k < npixels; k++)
        {
            out[k * ch + 0] = work[k * 3 + 0];
            out[k * ch + 1] = work[k * 3 + 1];
            out[k * ch + 2] = work[k * 3 + 2];
            out[k * ch + 3] = in[k * ch + 3];  /* preserve alpha */
        }
    }
    else
    {
        /* Standard density to transmittance (T = 10^(-D))
         * This creates a "negative-like" output where:
         * - High density (bright scene) -> low transmittance (dark in negative)
         * - Low density (dark scene) -> high transmittance (bright in negative)
         * Use negadoctor or invert module to get positive image.
         */
        const float LN10 = 2.302585093f;
        DT_OMP_FOR()
        for(size_t k = 0; k < npixels; k++)
        {
            float D_r = work[k * 3 + 0];
            float D_g = work[k * 3 + 1];
            float D_b = work[k * 3 + 2];
            
            /* T = 10^(-D) = exp(-D * ln(10)) */
            out[k * ch + 0] = expf(-D_r * LN10);
            out[k * ch + 1] = expf(-D_g * LN10);
            out[k * ch + 2] = expf(-D_b * LN10);
            out[k * ch + 3] = in[k * ch + 3];  /* preserve alpha */
        }
    }

    dt_free_align(work);
}

/* Load RGB-to-spectral model if not already loaded */
static void ensure_rgb2spec_loaded(void)
{
    if(g_rgb2spec_loaded) return;
    
    /* Try to load from user config directory */
    const char *config_dir = g_get_user_config_dir();
    if(config_dir)
    {
        char path[512];
        snprintf(path, sizeof(path), "%s/darktable/filmsim/jakob-and-hanika-2019-srgb.coeff", config_dir);
        if(filmsim_rgb2spec_load(&g_rgb2spec_model, path))
        {
            g_rgb2spec_loaded = TRUE;
            return;
        }
    }
    
    /* Try common system paths */
    const char *system_paths[] = {
        "/usr/share/cppfilmsim/data/jakob-and-hanika-2019-srgb.coeff",
        "/usr/local/share/cppfilmsim/data/jakob-and-hanika-2019-srgb.coeff",
        NULL
    };
    
    for(int i = 0; system_paths[i]; i++)
    {
        if(filmsim_rgb2spec_load(&g_rgb2spec_model, system_paths[i]))
        {
            g_rgb2spec_loaded = TRUE;
            return;
        }
    }
    
    fprintf(stderr, "[filmsim] RGB-to-spectral coefficients not found. Spectral mode will use fallback.\n");
}

void commit_params(struct dt_iop_module_t *self,
                   dt_iop_params_t *p1,
                   dt_dev_pixelpipe_t *pipe,
                   dt_dev_pixelpipe_iop_t *piece)
{
    dt_iop_filmsim_params_t *p = (dt_iop_filmsim_params_t *)p1;
    dt_iop_filmsim_data_t *data = (dt_iop_filmsim_data_t *)piece->data;

    data->push_pull = p->push_pull;
    data->enable_halation = p->enable_halation;
    data->enable_dir = p->enable_dir;
    data->enable_dye_mix = p->enable_dye_mix;
    data->enable_spectral = p->enable_spectral;
    data->halation_strength = p->halation_strength;
    data->dir_strength = p->dir_strength;
    data->pixel_size = p->pixel_size;

    int film_idx = CLAMP(p->film_profile, 0, FILMSIM_NUM_FILMS - 1);
    const char *film_name = filmsim_film_names[film_idx];
    const char *json_data = filmsim_get_embedded_profile(film_name);
    
    if(json_data)
    {
        filmsim_load_profile_from_json(&data->film, json_data);
    }

    /* Check if film has spectral data */
    data->has_spectral_data = filmsim_has_spectral_data(&data->film);
    
    /* Build spectral LUTs if spectral mode is enabled and data is available */
    if(data->enable_spectral && data->has_spectral_data)
    {
        ensure_rgb2spec_loaded();
        filmsim_spectral_free_luts(&data->spectral_lut);
        filmsim_spectral_build_luts(&data->spectral_lut, &data->film, 
                                     g_rgb2spec_loaded ? &g_rgb2spec_model : NULL);
    }

    build_film_lut(data);
}

void init_pipe(struct dt_iop_module_t *self,
               dt_dev_pixelpipe_t *pipe,
               dt_dev_pixelpipe_iop_t *piece)
{
    piece->data = calloc(1, sizeof(dt_iop_filmsim_data_t));
}

void cleanup_pipe(struct dt_iop_module_t *self,
                  dt_dev_pixelpipe_t *pipe,
                  dt_dev_pixelpipe_iop_t *piece)
{
    dt_iop_filmsim_data_t *data = (dt_iop_filmsim_data_t *)piece->data;
    if(data)
    {
        filmsim_spectral_free_luts(&data->spectral_lut);
        free(data);
    }
    piece->data = NULL;
}

static void film_callback(GtkWidget *widget, dt_iop_module_t *self)
{
    if(darktable.gui->reset) return;
    dt_iop_filmsim_params_t *p = (dt_iop_filmsim_params_t *)self->params;
    p->film_profile = dt_bauhaus_combobox_get(widget);
    dt_dev_add_history_item(darktable.develop, self, TRUE);
}

void gui_init(struct dt_iop_module_t *self)
{
    dt_iop_filmsim_gui_data_t *g = IOP_GUI_ALLOC(filmsim);

    /* Film profile combobox */
    g->film_combo = dt_bauhaus_combobox_new(self);
    dt_bauhaus_widget_set_label(g->film_combo, NULL, N_("film profile"));
    
    dt_bauhaus_combobox_add(g->film_combo, _("Kodak Portra 400"));
    dt_bauhaus_combobox_add(g->film_combo, _("Kodak Ektar 100"));
    dt_bauhaus_combobox_add(g->film_combo, _("Kodak Ultramax 400"));
    dt_bauhaus_combobox_add(g->film_combo, _("Kodak Vision3 500T"));
    dt_bauhaus_combobox_add(g->film_combo, _("Kodak Kodachrome 64"));
    dt_bauhaus_combobox_add(g->film_combo, _("Kodak Ektachrome P 1600"));
    dt_bauhaus_combobox_add(g->film_combo, _("Kodak Royal Gold 1000"));
    dt_bauhaus_combobox_add(g->film_combo, _("Kodak Portra 100T"));
    dt_bauhaus_combobox_add(g->film_combo, _("Kodak T-Max 400"));
    dt_bauhaus_combobox_add(g->film_combo, _("Fuji Pro H 400"));
    dt_bauhaus_combobox_add(g->film_combo, _("Fuji Superia 1600"));
    dt_bauhaus_combobox_add(g->film_combo, _("Fuji Superia Reala 100"));
    dt_bauhaus_combobox_add(g->film_combo, _("Fuji Fujichrome Velvia 100"));
    dt_bauhaus_combobox_add(g->film_combo, _("Fuji Fujichrome Astia F 100"));
    dt_bauhaus_combobox_add(g->film_combo, _("Fuji Neopan Acros II 100"));
    dt_bauhaus_combobox_add(g->film_combo, _("Cinestill 800T"));
    dt_bauhaus_combobox_add(g->film_combo, _("Agfacolor Optima II 100"));
    dt_bauhaus_combobox_add(g->film_combo, _("Agfacolor Portrait XPS 160"));
    
    gtk_widget_set_tooltip_text(g->film_combo, _("select analog film stock to simulate"));
    g_signal_connect(G_OBJECT(g->film_combo), "value-changed", G_CALLBACK(film_callback), self);

    /* Create widget container with combobox first */
    self->widget = dt_gui_vbox(g->film_combo);

    /* Push/pull slider */
    g->push_pull_slider = dt_bauhaus_slider_from_params(self, "push_pull");
    dt_bauhaus_slider_set_format(g->push_pull_slider, _("%.1f stops"));
    gtk_widget_set_tooltip_text(g->push_pull_slider, _("simulate push/pull development"));

    /* Halation toggle and strength */
    g->halation_toggle = dt_bauhaus_toggle_from_params(self, "enable_halation");
    gtk_widget_set_tooltip_text(g->halation_toggle, _("enable halation (light scatter from film base)"));

    g->halation_strength_slider = dt_bauhaus_slider_from_params(self, "halation_strength");
    gtk_widget_set_tooltip_text(g->halation_strength_slider, _("strength of halation effect"));

    /* DIR toggle and strength */
    g->dir_toggle = dt_bauhaus_toggle_from_params(self, "enable_dir");
    gtk_widget_set_tooltip_text(g->dir_toggle, _("enable DIR coupler edge enhancement"));

    g->dir_strength_slider = dt_bauhaus_slider_from_params(self, "dir_strength");
    gtk_widget_set_tooltip_text(g->dir_strength_slider, _("strength of DIR edge enhancement"));

    /* Dye mixing toggle */
    g->dye_mix_toggle = dt_bauhaus_toggle_from_params(self, "enable_dye_mix");
    gtk_widget_set_tooltip_text(g->dye_mix_toggle, _("enable cross-layer dye contamination"));

    /* Spectral mode toggle */
    g->spectral_toggle = dt_bauhaus_toggle_from_params(self, "enable_spectral");
    gtk_widget_set_tooltip_text(g->spectral_toggle, _("use spectral simulation for more accurate color (requires spectral data in film profile)"));

    /* Pixel size slider */
    g->pixel_size_slider = dt_bauhaus_slider_from_params(self, "pixel_size");
    dt_bauhaus_slider_set_format(g->pixel_size_slider, _("%.1f µm"));
    gtk_widget_set_tooltip_text(g->pixel_size_slider, _("pixel size in micrometers (affects halation/DIR spread)"));
}

void gui_update(struct dt_iop_module_t *self)
{
    dt_iop_filmsim_params_t *p = (dt_iop_filmsim_params_t *)self->params;
    dt_iop_filmsim_gui_data_t *g = (dt_iop_filmsim_gui_data_t *)self->gui_data;

    dt_bauhaus_combobox_set(g->film_combo, p->film_profile);
    /* sliders and toggles auto-update from from_params */
}

void gui_cleanup(struct dt_iop_module_t *self)
{
    free(self->gui_data);
    self->gui_data = NULL;
}
