#!/bin/bash
#
# Build darktable with the film simulation plugin
#
# This script clones darktable, integrates the filmsim plugin,
# builds everything, and creates a ready-to-run installation.
#
# Usage:
#   ./scripts/build_darktable_with_filmsim.sh [options]
#
# Output:
#   build/darktable_custom/ - Complete runnable darktable build
#

set -e

# ============================================================================
# Configuration
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="${PROJECT_ROOT}/build"
DARKTABLE_SRC="${BUILD_DIR}/darktable"
DARKTABLE_BUILD="${BUILD_DIR}/darktable_build"
OUTPUT_DIR="${BUILD_DIR}/darktable_custom"
FILM_PROFILES_DIR="$HOME/.config/darktable/filmsim"

# Darktable version/branch to use
DARKTABLE_REPO="https://github.com/darktable-org/darktable.git"
DARKTABLE_BRANCH="master"

# Build settings
NUM_JOBS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
CLEAN_BUILD=false
UPGRADE_DARKTABLE=false

# ============================================================================
# Parse arguments
# ============================================================================

while [[ $# -gt 0 ]]; do
    case "$1" in
        --clean)
            CLEAN_BUILD=true
            shift
            ;;
        --upgrade-darktable)
            UPGRADE_DARKTABLE=true
            shift
            ;;
        --jobs)
            NUM_JOBS="$2"
            shift 2
            ;;
        -j*)
            NUM_JOBS="${1#-j}"
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [options]"
            echo ""
            echo "Options:"
            echo "  --clean              Remove existing build and start fresh"
            echo "  --upgrade-darktable  Pull latest darktable from upstream (risky)"
            echo "  --jobs N             Number of parallel build jobs (default: auto)"
            echo ""
            echo "Output: ${OUTPUT_DIR}/"
            echo ""
            echo "The script is idempotent - safe to run multiple times."
            echo "It will reuse existing darktable source unless --upgrade-darktable is used."
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# ============================================================================
# Helper functions
# ============================================================================

log() {
    echo "[$(date '+%H:%M:%S')] $*"
}

error() {
    echo "[ERROR] $*" >&2
    exit 1
}

check_dependencies() {
    log "Checking build dependencies..."
    
    local missing=()
    
    # Required tools
    for cmd in git cmake gcc pkg-config; do
        if ! command -v "$cmd" &>/dev/null; then
            missing+=("$cmd")
        fi
    done
    
    # Required libraries (check via pkg-config)
    for lib in gtk+-3.0 glib-2.0 json-glib-1.0 fftw3f; do
        if ! pkg-config --exists "$lib" 2>/dev/null; then
            missing+=("$lib")
        fi
    done
    
    if [[ ${#missing[@]} -gt 0 ]]; then
        error "Missing dependencies: ${missing[*]}

Install on Ubuntu/Debian:
  sudo apt install git cmake gcc g++ pkg-config \\
    libgtk-3-dev libglib2.0-dev libjson-glib-dev libfftw3-dev \\
    libxml2-dev libsqlite3-dev libcurl4-openssl-dev \\
    libpng-dev libjpeg-dev libtiff-dev libexiv2-dev \\
    liblensfun-dev libpugixml-dev libgphoto2-dev \\
    liblcms2-dev libcolord-dev libcolord-gtk-dev \\
    libgraphicsmagick1-dev libopenexr-dev \\
    intltool xsltproc

Install on Fedora:
  sudo dnf install git cmake gcc gcc-c++ pkgconfig \\
    gtk3-devel glib2-devel json-glib-devel fftw-devel \\
    libxml2-devel sqlite-devel libcurl-devel \\
    libpng-devel libjpeg-turbo-devel libtiff-devel exiv2-devel \\
    lensfun-devel pugixml-devel libgphoto2-devel \\
    lcms2-devel colord-devel colord-gtk-devel \\
    GraphicsMagick-devel openexr-devel \\
    intltool libxslt

Install on macOS (Homebrew):
  brew install cmake pkg-config gtk+3 json-glib fftw \\
    libxml2 sqlite curl libpng jpeg libtiff exiv2 \\
    lensfun pugixml libgphoto2 little-cms2 colord \\
    graphicsmagick openexr intltool libxslt libomp"
    fi
    
    log "All dependencies found"
}

# ============================================================================
# Main build process
# ============================================================================

main() {
    log "=================================================="
    log "Building darktable with film simulation plugin"
    log "=================================================="
    log "Project root: ${PROJECT_ROOT}"
    log "Build directory: ${BUILD_DIR}"
    log "Output directory: ${OUTPUT_DIR}"
    log "Parallel jobs: ${NUM_JOBS}"
    log ""
    
    check_dependencies
    
    # Clean if requested
    if [[ "$CLEAN_BUILD" == "true" ]]; then
        log "Cleaning previous build..."
        rm -rf "$BUILD_DIR"
    fi
    
    mkdir -p "$BUILD_DIR"
    
    # ========================================================================
    # Step 1: Clone darktable (only if not present, or upgrade requested)
    # ========================================================================
    
    if [[ -d "$DARKTABLE_SRC/.git" ]]; then
        if [[ "$UPGRADE_DARKTABLE" == "true" ]]; then
            log "Upgrading darktable source (--upgrade-darktable)..."
            cd "$DARKTABLE_SRC"
            git fetch origin
            git checkout "$DARKTABLE_BRANCH"
            git pull origin "$DARKTABLE_BRANCH"
            git submodule update --init --recursive
            log "Darktable upgraded to latest"
        else
            log "Using existing darktable source (use --upgrade-darktable to update)"
            # Ensure submodules are initialized even for existing source
            cd "$DARKTABLE_SRC"
            if [[ ! -f "src/external/rawspeed/CMakeLists.txt" ]]; then
                log "Initializing submodules..."
                git submodule update --init --recursive
            fi
        fi
    else
        log "Cloning darktable repository (this may take a few minutes)..."
        rm -rf "$DARKTABLE_SRC"
        git clone --recursive --depth 1 --branch "$DARKTABLE_BRANCH" "$DARKTABLE_REPO" "$DARKTABLE_SRC"
    fi
    
    # ========================================================================
    # Step 2: Copy plugin files to darktable source
    # ========================================================================
    
    log "Copying filmsim plugin files..."
    
    local IOP_DIR="${DARKTABLE_SRC}/src/iop"
    
    cp "${PROJECT_ROOT}/src/filmsim.c" "$IOP_DIR/"
    cp "${PROJECT_ROOT}/src/filmsim_common.h" "$IOP_DIR/"
    cp "${PROJECT_ROOT}/src/filmsim_curves.c" "$IOP_DIR/"
    cp "${PROJECT_ROOT}/src/filmsim_halation.c" "$IOP_DIR/"
    cp "${PROJECT_ROOT}/src/filmsim_dir.c" "$IOP_DIR/"
    cp "${PROJECT_ROOT}/src/filmsim_curves_app.c" "$IOP_DIR/"
    cp "${PROJECT_ROOT}/src/filmsim_data.c" "$IOP_DIR/"
    cp "${PROJECT_ROOT}/src/filmsim_spectral.c" "$IOP_DIR/"
    cp "${PROJECT_ROOT}/src/filmsim_spectral.h" "$IOP_DIR/"
    
    log "Plugin files copied to ${IOP_DIR}/"
    
    # ========================================================================
    # Step 3: Patch CMakeLists.txt to add filmsim module (idempotent)
    # ========================================================================
    
    log "Patching CMakeLists.txt..."
    
    local CMAKE_FILE="${IOP_DIR}/CMakeLists.txt"
    
    # Check if already patched
    if grep -q "add_iop(filmsim" "$CMAKE_FILE" 2>/dev/null; then
        log "CMakeLists.txt already patched, skipping..."
    else
        # Find the line with grain module and add filmsim after it
        local PATCH_LINE='add_iop(filmsim "filmsim.c" "filmsim_curves.c" "filmsim_halation.c" "filmsim_dir.c" "filmsim_curves_app.c" "filmsim_data.c" "filmsim_spectral.c")'
        # Note: must use plain signature (no PRIVATE/PUBLIC) to match add_iop macro style
        local LINK_LINE='target_link_libraries(filmsim fftw3f)'
        
        # Find line number of grain module
        local grain_line=$(grep -n 'add_iop(grain' "$CMAKE_FILE" | head -1 | cut -d: -f1)
        
        if [[ -z "$grain_line" ]]; then
            # If no grain module found, add at end
            log "Warning: grain module not found, adding filmsim at end of file"
            echo "" >> "$CMAKE_FILE"
            echo "# Film Simulation plugin" >> "$CMAKE_FILE"
            echo "$PATCH_LINE" >> "$CMAKE_FILE"
            echo "$LINK_LINE" >> "$CMAKE_FILE"
        else
            # Insert after grain module line
            local temp_file=$(mktemp)
            head -n "$grain_line" "$CMAKE_FILE" > "$temp_file"
            echo "" >> "$temp_file"
            echo "# Film Simulation plugin" >> "$temp_file"
            echo "$PATCH_LINE" >> "$temp_file"
            echo "$LINK_LINE" >> "$temp_file"
            tail -n +$((grain_line + 1)) "$CMAKE_FILE" >> "$temp_file"
            mv "$temp_file" "$CMAKE_FILE"
        fi
        
        log "CMakeLists.txt patched"
    fi
    
    # ========================================================================
    # Step 4: Patch iop_order.c to add filmsim to pipeline (idempotent)
    # ========================================================================
    
    log "Patching iop_order.c..."
    
    local IOP_ORDER_FILE="${DARKTABLE_SRC}/src/common/iop_order.c"
    
    if grep -q '"filmsim"' "$IOP_ORDER_FILE" 2>/dev/null; then
        log "iop_order.c already patched, skipping..."
    else
        # There are multiple iop_order tables in this file (legacy_order, v30_order, 
        # v50_order, v30_jpg_order, v50_jpg_order). We need to add filmsim to ALL of them.
        # Strategy: Find each occurrence of "grain" and insert filmsim after it.
        # The grain module appears in each table, so this adds filmsim to all tables.
        # We use a unique order number that fits in the creative module range.
        
        # Use sed to insert filmsim after every "grain" line
        # Handle both formats: { {56.0f } and { { 65.0f } (with extra space)
        # The order value 65.5 places it between grain (65.0 in v30/v50) and soften (66.0)
        sed -i 's/\({ *{ *[0-9.]*f *}, "grain", 0 *},\?\)/\1\n  { {65.5f }, "filmsim", 0},/' "$IOP_ORDER_FILE"
        
        local count=$(grep -c '"filmsim"' "$IOP_ORDER_FILE")
        log "iop_order.c patched ($count entries added)"
    fi
    
    # ========================================================================
    # Step 5: Configure and build darktable
    # ========================================================================
    
    log "Configuring darktable build..."
    
    mkdir -p "$DARKTABLE_BUILD"
    cd "$DARKTABLE_BUILD"
    
    # Detect OS for platform-specific options
    local CMAKE_EXTRA_OPTS=""
    if [[ "$(uname)" == "Darwin" ]]; then
        # macOS: help find OpenMP from Homebrew
        if [[ -d "$(brew --prefix)/opt/libomp" ]]; then
            export OpenMP_ROOT="$(brew --prefix)/opt/libomp"
        fi
    fi
    
    cmake "$DARKTABLE_SRC" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="$OUTPUT_DIR" \
        -DBUILD_USERMANUAL=OFF \
        -DBUILD_TESTING=OFF \
        -DUSE_OPENCL=ON \
        -DUSE_LUA=ON \
        $CMAKE_EXTRA_OPTS
    
    log "Building darktable (this may take a while)..."
    
    cmake --build . -j "$NUM_JOBS"
    
    # ========================================================================
    # Step 6: Install to output directory
    # ========================================================================
    
    log "Installing to ${OUTPUT_DIR}..."
    
    cmake --install .
    
    # ========================================================================
    # Step 7: Install film profiles
    # ========================================================================
    
    log "Installing film profiles..."
    
    mkdir -p "$FILM_PROFILES_DIR"
    cp "${PROJECT_ROOT}/data/film_profiles/"*.json "$FILM_PROFILES_DIR/"
    
    # Copy RGB-to-spectral coefficients if present
    if [[ -f "${PROJECT_ROOT}/data/jakob-and-hanika-2019-srgb.coeff" ]]; then
        cp "${PROJECT_ROOT}/data/jakob-and-hanika-2019-srgb.coeff" "$FILM_PROFILES_DIR/"
        log "Installed RGB-to-spectral coefficients"
    fi
    
    local profile_count=$(ls -1 "${FILM_PROFILES_DIR}"/*.json 2>/dev/null | wc -l)
    log "Installed ${profile_count} film profiles to ${FILM_PROFILES_DIR}/"
    
    # ========================================================================
    # Done!
    # ========================================================================
    
    log ""
    log "=================================================="
    log "Build complete!"
    log "=================================================="
    log ""
    log "To run darktable with film simulation:"
    log "  ${OUTPUT_DIR}/bin/darktable"
    log ""
    log "The film simulation module appears in the 'effects' group"
    log "in the darkroom. It outputs a negative image that works"
    log "well with darktable's negadoctor module."
    log ""
    log "Available film profiles (${profile_count}):"
    ls -1 "${FILM_PROFILES_DIR}"/*.json 2>/dev/null | xargs -I{} basename {} .json | sed 's/^/  - /'
    log ""
}

# Run main
main "$@"
