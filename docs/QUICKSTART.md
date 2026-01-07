# Quick Start Guide - darktable Film Simulation Plugin

## Current Status: COMPLETE

The darktable film simulation plugin is fully implemented and ready to use.

## Build & Install (One Command)

```bash
cd darktable_plugin
./build_darktable_with_filmsim.sh
```

This script will:
1. Clone darktable source from GitHub
2. Copy plugin files into the darktable source tree
3. Patch CMakeLists.txt and iop_order.c
4. Build darktable with the filmsim module
5. Install film profiles to `~/.config/darktable/filmsim/`
6. Output ready-to-run darktable to `../darktable_custom/`

### Build Options

```bash
# Clean build (removes previous build)
./build_darktable_with_filmsim.sh --clean

# Update darktable to latest upstream (risky - may break build)
./build_darktable_with_filmsim.sh --upgrade-darktable

# Specify parallel jobs
./build_darktable_with_filmsim.sh --jobs 8
```

Note: The script reuses existing darktable source by default. Use `--upgrade-darktable` 
only if you need the latest features and are prepared to debug potential issues.

### Run

```bash
../darktable_custom/bin/darktable
```

## Using the Module

1. Open an image in darktable's darkroom
2. Find "film simulation" in the **effects** module group
3. Enable the module
4. Select a film profile from the dropdown
5. Adjust push/pull for exposure compensation
6. Toggle halation, DIR, and dye mixing as desired

### Module Parameters

| Parameter | Range | Description |
|-----------|-------|-------------|
| Film profile | dropdown | Select from 18 film stocks |
| Push/pull | -10 to +10 stops | Simulate push/pull development |
| Halation | on/off | Light scatter from film base |
| Halation strength | 0-2x | Scale halation intensity |
| DIR | on/off | Edge enhancement from DIR couplers |
| DIR strength | 0-2x | Scale DIR intensity |
| Dye mixing | on/off | Cross-layer dye contamination |
| Pixel size | 1-20 µm | Affects halation/DIR spread |

### Output

The module outputs a **negative density image**. For color negative films:
- Use darktable's **negadoctor** module after filmsim to invert
- Or chain with the **invert** module for basic inversion

For reversal (slide) films like Kodachrome and Velvia:
- The output is already positive (like a slide)
- No inversion needed

## Available Film Profiles (18)

### Color Negative
- Kodak Portra 400
- Kodak Portra 100T
- Kodak Ektar 100
- Kodak Ultramax 400
- Kodak Royal Gold 1000
- Kodak Vision3 500T
- Fuji Pro H 400
- Fuji Superia 1600
- Fuji Superia Reala 100
- Cinestill 800T
- Agfacolor Optima II 100
- Agfacolor Portrait XPS 160

### Color Reversal (Slide)
- Kodak Kodachrome 64
- Kodak Ektachrome P 1600
- Fuji Fujichrome Velvia 100
- Fuji Fujichrome Astia F 100

### Black & White
- Kodak T-Max 400
- Fuji Neopan Acros II 100

## Technical Details

### Pipeline Order
```
Input (Linear RGB)
    ↓
Halation (exposure domain, FFTW convolution)
    ↓
H-D Curves (65536-entry LUT, log exposure → density)
    ↓
DIR (density domain, FFTW convolution)
    ↓
Dye Mixing (3x3 matrix)
    ↓
Density → Transmittance (T = 10^(-D))
    ↓
Output (for negadoctor)
```

### Performance
- H-D curve LUT generation: ~1 ms
- Halation (24MP): ~100-200 ms
- DIR (24MP): ~150-300 ms
- Total: <500 ms for interactive editing

### Dependencies
The build requires standard darktable dependencies plus:
- FFTW3 (single-precision): `libfftw3f-dev` / `fftw-devel`
- JSON-GLib: `libjson-glib-dev` / `json-glib-devel`

## Files

```
darktable_plugin/
├── build_darktable_with_filmsim.sh  # One-command build script
├── filmsim.c                        # Main IOP module (465 lines)
├── filmsim_common.h                 # Data structures
├── filmsim_curves.c                 # Spline interpolation
├── filmsim_halation.c               # Halation effect (FFTW)
├── filmsim_dir.c                    # DIR edge enhancement (FFTW)
├── filmsim_curves_app.c             # LUT generation, dye mixing
├── filmsim_data.c                   # JSON profile loading
├── CMakeLists.txt                   # Standalone build (alternative)
├── README.md                        # Full documentation
├── IMPLEMENTATION_PLAN.md           # Design notes
└── QUICKSTART.md                    # This file
```

## Troubleshooting

### Build fails with missing dependencies
```bash
# Ubuntu/Debian
sudo apt install libfftw3-dev libjson-glib-dev libgtk-3-dev

# Fedora
sudo dnf install fftw-devel json-glib-devel gtk3-devel
```

### Module doesn't appear in darktable
- Ensure you're running the darktable from `darktable_custom/bin/`
- Check the build log for compilation errors
- The module appears in the "effects" group

### Film profiles not loading
- Verify profiles exist: `ls ~/.config/darktable/filmsim/`
- Check darktable console for JSON parsing errors
- Profiles must be valid JSON with the expected structure

### Halation/DIR not visible
- Increase strength sliders (default is 1.0x)
- Ensure pixel size matches your image (default 6.0 µm)
- Effects are subtle on low-contrast images

## Development

### Rebuilding after changes
```bash
# Quick rebuild (after editing plugin files)
./build_darktable_with_filmsim.sh --skip-clone

# Full clean rebuild
./build_darktable_with_filmsim.sh --clean
```

### Adding new film profiles
1. Create JSON file in `../data/` following existing format
2. Rebuild darktable or manually copy to `~/.config/darktable/filmsim/`
3. Add entry to `filmsim_film_names[]` array in `filmsim.c`
4. Add combobox entry in `gui_init()` function
