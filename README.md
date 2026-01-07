# darktable Film Simulation Plugin

A native darktable plugin that simulates analog film characteristics with physically-based models. This plugin implements the complete film photography pipeline including H-D curves, halation, DIR (Development Inhibitor Releasing) couplers, dye mixing, and optional spectral simulation.

## Features

- **18 Film Profiles**: Kodak, Fuji, Agfa, and Cinestill films
- **H-D Curves**: Hurter-Driffield characteristic curves with push/pull development
- **Halation**: Light scatter from film base (wavelength-dependent)
- **DIR Couplers**: Edge enhancement from developer inhibitor diffusion
- **Dye Mixing**: Cross-layer dye contamination simulation
- **Spectral Mode**: Optional physically-accurate spectral simulation using Jakob & Hanika 2019 RGB-to-spectral conversion

## Supported Film Stocks

### Color Negative
- Kodak Portra 400, Portra 100T, Ektar 100, Ultramax 400, Vision3 500T, Royal Gold 1000
- Fuji Pro H 400, Superia 1600, Superia Reala 100
- Cinestill 800T
- Agfacolor Optima II 100, Portrait XPS 160

### Color Reversal (Slide)
- Kodak Kodachrome 64, Ektachrome P 1600
- Fuji Fujichrome Velvia 100, Astia F 100

### Black & White
- Kodak T-Max 400
- Fuji Neopan Acros II 100

## Installation

### Pre-built Releases (Recommended)

Download the latest release from the [Releases page](../../releases).

#### macOS
1. Download `darktable-filmsim-macos-arm64.dmg`
2. Open the DMG and drag **darktable (Film Simulation)** to your Applications folder
3. Film profiles and spectral data are bundled inside the app - no additional setup needed!

#### Linux (AppImage)
1. Download `darktable-filmsim-linux-x86_64.AppImage`
2. Make it executable: `chmod +x darktable-filmsim-linux-x86_64.AppImage`
3. Run it: `./darktable-filmsim-linux-x86_64.AppImage`
4. Film profiles are bundled - no additional setup needed!

#### Linux (Tarball)
1. Download and extract `darktable-filmsim-linux-x86_64.tar.gz`
2. **Important**: Copy film profiles to your darktable config:
   ```bash
   mkdir -p ~/.config/darktable/filmsim
   cp darktable_custom/share/darktable/filmsim/* ~/.config/darktable/filmsim/
   ```
3. Run: `./darktable_custom/bin/darktable`

### Building from Source

If you prefer to build from source, you'll need development dependencies.

#### Prerequisites

<details>
<summary>Ubuntu/Debian</summary>

```bash
sudo apt install git cmake gcc g++ pkg-config \
  libgtk-3-dev libglib2.0-dev libjson-glib-dev libfftw3-dev \
  libxml2-dev libsqlite3-dev libcurl4-openssl-dev \
  libpng-dev libjpeg-dev libtiff-dev libexiv2-dev \
  liblensfun-dev libpugixml-dev libgphoto2-dev \
  liblcms2-dev libcolord-dev libcolord-gtk-dev \
  libgraphicsmagick1-dev libopenexr-dev \
  librsvg2-dev libwebp-dev \
  intltool xsltproc
```
</details>

<details>
<summary>Fedora</summary>

```bash
sudo dnf install git cmake gcc gcc-c++ pkgconfig \
  gtk3-devel glib2-devel json-glib-devel fftw-devel \
  libxml2-devel sqlite-devel libcurl-devel \
  libpng-devel libjpeg-turbo-devel libtiff-devel exiv2-devel \
  lensfun-devel pugixml-devel libgphoto2-devel \
  lcms2-devel colord-devel colord-gtk-devel \
  GraphicsMagick-devel openexr-devel \
  librsvg2-devel libwebp-devel \
  intltool libxslt
```
</details>

<details>
<summary>macOS (Homebrew)</summary>

```bash
brew install cmake pkg-config gtk+3 json-glib fftw \
  libxml2 sqlite curl libpng jpeg-turbo libtiff exiv2 \
  lensfun pugixml libgphoto2 little-cms2 \
  graphicsmagick openexr intltool libxslt libomp \
  librsvg webp
```
</details>

#### Build Steps

```bash
# Clone the repository (with Git LFS for large files)
git lfs install
git clone https://github.com/YOUR_USERNAME/darktable_filmsim.git
cd darktable_filmsim

# Build darktable with the film simulation plugin
./scripts/build_darktable_with_filmsim.sh

# Install film profiles (REQUIRED for source builds)
mkdir -p ~/.config/darktable/filmsim
cp data/film_profiles/*.json ~/.config/darktable/filmsim/
cp data/jakob-and-hanika-2019-srgb.coeff ~/.config/darktable/filmsim/

# Run darktable
./build/darktable_custom/bin/darktable
```

Build options:
- `--clean`: Remove existing build and start fresh
- `--upgrade-darktable`: Pull latest darktable from upstream
- `--jobs N` or `-jN`: Set number of parallel build jobs

## Usage

1. Open an image in darktable's darkroom
2. Find "film simulation" in the **effects** module group
3. Select a film profile from the dropdown
4. Adjust parameters:
   - **Push/Pull**: Simulate over/under development (-10 to +10 stops)
   - **Halation**: Enable/disable and adjust strength
   - **DIR**: Enable/disable and adjust edge enhancement strength
   - **Dye Mixing**: Enable/disable cross-layer dye effects
   - **Spectral Mode**: Enable for more accurate color (slower)
   - **Pixel Size**: Physical pixel size in micrometers (affects halation/DIR spread)

5. The module outputs a **film negative** - use darktable's **negadoctor** module afterwards to invert and color-correct

### Recommended Workflow

1. Enable film simulation module
2. Select your desired film stock
3. Adjust push/pull if desired
4. Add negadoctor module after film simulation
5. Use negadoctor's color picker on the film border (if present) to set D-min
6. Adjust negadoctor settings for final look

## Film Profile Location

The plugin looks for film profiles in this order:
1. **Bundled location** (for AppImage/DMG): `<app>/share/darktable/filmsim/`
2. **User config**: `~/.config/darktable/filmsim/`

Each film profile is a JSON file containing:
- H-D characteristic curves
- Spectral sensitivity data
- Halation parameters
- DIR coupler settings
- Dye mixing matrix

The `jakob-and-hanika-2019-srgb.coeff` file (~9MB) contains RGB-to-spectral conversion coefficients for the optional spectral mode.

## Technical Details

### Pipeline Order

```
Input (Linear RGB)
    ↓
[Spectral Exposure] (optional)
    ↓
Halation (FFTW convolution)
    ↓
H-D Curves (65536-entry LUT)
    ↓
DIR (FFTW convolution)
    ↓
Dye Mixing (3x3 matrix)
    ↓
[Spectral Dye Model] or Density→Transmittance
    ↓
Output (Negative for negadoctor)
```

### Spectral Mode

When enabled, spectral mode uses:
- Jakob & Hanika 2019 RGB-to-spectral polynomial model
- 32³ LUTs for fast spectral integration
- Film spectral sensitivity curves
- Wavelength-dependent dye absorption (kappa functions)
- 5000K scanner illuminant simulation

This provides more physically accurate color reproduction but is slower than the standard mode.

## File Structure

```
darktable_filmsim/
├── src/                    # Plugin source code
│   ├── filmsim.c          # Main darktable IOP module
│   ├── filmsim_common.h   # Data structures
│   ├── filmsim_curves.c   # Spline interpolation
│   ├── filmsim_data.c     # JSON profile loading
│   ├── filmsim_halation.c # Halation effect (FFTW)
│   ├── filmsim_dir.c      # DIR effect (FFTW)
│   ├── filmsim_spectral.c # Spectral processing
│   └── filmsim_spectral.h
├── data/
│   ├── film_profiles/     # 18 JSON film profiles
│   └── jakob-and-hanika-2019-srgb.coeff  # Spectral coefficients (LFS)
├── scripts/
│   └── build_darktable_with_filmsim.sh
└── docs/
    ├── IMPLEMENTATION_PLAN.md
    ├── QUICKSTART.md
    └── STATUS.md
```

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

This is the same license as darktable, ensuring compatibility.

## Acknowledgments

- [darktable](https://www.darktable.org/) - The open source photography workflow application
- Jakob & Hanika 2019 - "A Low-Dimensional Function Space for Efficient Spectral Upsampling"
- Film characteristic data derived from manufacturer specifications and measurements

## Contributing

Contributions are welcome! Please feel free to submit issues and pull requests.

## Related Projects

- [cppfilmsim](https://github.com/YOUR_USERNAME/cppfilmsim) - Standalone C++ film simulation tool (this plugin's origin)
