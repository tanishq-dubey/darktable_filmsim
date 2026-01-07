/*
 * Film Simulation Plugin for darktable
 * Common data structures and definitions
 */

#pragma once

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Maximum number of curve points per film profile */
#define FILMSIM_MAX_CURVE_POINTS 256

/* Maximum number of spectral bands */
#define FILMSIM_MAX_SPECTRAL_BANDS 128

/* Film types */
typedef enum {
    FILMSIM_FILM_TYPE_COLOR_NEGATIVE = 0,
    FILMSIM_FILM_TYPE_COLOR_POSITIVE = 1,  /* Reversal/Slide */
    FILMSIM_FILM_TYPE_MONOCHROME_NEGATIVE = 2
} filmsim_film_type_t;

/* Film profile information */
typedef struct {
    char name[128];
    char description[256];
    int format_mm;
    filmsim_film_type_t film_type;
} filmsim_film_info_t;

/* Calibration data */
typedef struct {
    int iso;
    double middle_gray_logE;
} filmsim_calibration_t;

/* Halation parameters */
typedef struct {
    struct {
        double r, g, b;
    } strength;
    struct {
        double r, g, b;
    } size_um;
} filmsim_halation_t;

/* Coupler parameters (DIR and dye mixing) */
typedef struct {
    double saturation_amount;
    double dir_amount_rgb[3];  /* C, M, Y */
    double dir_diffusion_um;
    double dir_diffusion_interlayer;
    
    /* Dye mixing matrix [3][3] */
    /* Row = output dye, Col = input dye */
    /* [0] = C output from C,M,Y inputs */
    /* [1] = M output from C,M,Y inputs */
    /* [2] = Y output from C,M,Y inputs */
    double dye_mix_matrix[3][3];
} filmsim_couplers_t;

/* Curve point (color film) */
typedef struct {
    double d;  /* log exposure */
    double r;  /* Red channel density/response */
    double g;  /* Green channel density/response */
    double b;  /* Blue channel density/response */
} filmsim_curve_point_rgb_t;

/* Curve point (monochrome film) */
typedef struct {
    double d;       /* log exposure */
    double density; /* Density */
} filmsim_curve_point_mono_t;

/* Spectral density point */
typedef struct {
    double wavelength;  /* nm */
    double density;     /* Optical density */
} filmsim_spectral_density_t;

/* Spectral sensitivity point (color) */
typedef struct {
    double wavelength;  /* nm */
    double y;  /* Yellow layer sensitivity */
    double m;  /* Magenta layer sensitivity */
    double c;  /* Cyan layer sensitivity */
} filmsim_spectral_sensitivity_rgb_t;

/* Spectral sensitivity point (mono) */
typedef struct {
    double wavelength;  /* nm */
    double sensitivity;
} filmsim_spectral_sensitivity_mono_t;

/* Complete film profile */
typedef struct {
    filmsim_film_info_t info;
    filmsim_calibration_t calibration;
    filmsim_halation_t halation;
    filmsim_couplers_t couplers;
    
    /* H-D curves (either RGB or mono, determined by film_type) */
    union {
        filmsim_curve_point_rgb_t rgb[FILMSIM_MAX_CURVE_POINTS];
        filmsim_curve_point_mono_t mono[FILMSIM_MAX_CURVE_POINTS];
    } curves;
    int num_curves;
    
    /* Spectral sensitivity (either RGB or mono) */
    union {
        filmsim_spectral_sensitivity_rgb_t rgb[FILMSIM_MAX_SPECTRAL_BANDS];
        filmsim_spectral_sensitivity_mono_t mono[FILMSIM_MAX_SPECTRAL_BANDS];
    } spectral_sensitivity;
    int num_spectral_sensitivity;
    
    /* Spectral density curves (negative films) */
    filmsim_spectral_density_t spectral_density_min[FILMSIM_MAX_SPECTRAL_BANDS];
    filmsim_spectral_density_t spectral_density_mid[FILMSIM_MAX_SPECTRAL_BANDS];
    int num_spectral_density_min;
    int num_spectral_density_mid;
    
    /* Per-dye spectral density (reversal/slide films) */
    filmsim_spectral_density_t spectral_density_y[FILMSIM_MAX_SPECTRAL_BANDS];
    filmsim_spectral_density_t spectral_density_m[FILMSIM_MAX_SPECTRAL_BANDS];
    filmsim_spectral_density_t spectral_density_c[FILMSIM_MAX_SPECTRAL_BANDS];
    int num_spectral_density_y;
    int num_spectral_density_m;
    int num_spectral_density_c;
} filmsim_film_profile_t;

/* Monotone cubic spline for curve interpolation */
typedef struct {
    double *x;       /* X coordinates (allocated) */
    double *y;       /* Y coordinates (allocated) */
    double *m;       /* Tangents (allocated) */
    int n;           /* Number of points */
} filmsim_spline_t;

/* Initialize spline (allocates memory) */
void filmsim_spline_init(filmsim_spline_t *spline);

/* Fit spline to points (must be sorted by x) */
void filmsim_spline_fit(filmsim_spline_t *spline, const double *x, const double *y, int n);

/* Interpolate value at x */
double filmsim_spline_interpolate(const filmsim_spline_t *spline, double x);

/* Free spline memory */
void filmsim_spline_free(filmsim_spline_t *spline);

/* Halation effect */
void filmsim_apply_halation(
    float *imageData,
    int width,
    int height,
    const filmsim_halation_t *halation,
    float pixelSizeUm
);

/* DIR effect */
void filmsim_apply_dir(
    float *imageData,
    int width,
    int height,
    const filmsim_couplers_t *couplers,
    float pixelSizeUm
);

/* Film profile application (H-D curves) */
void filmsim_apply_film_profile(
    float *imageData,
    int width,
    int height,
    const filmsim_film_profile_t *film,
    float pushPullStops
);

/* Dye mixing */
void filmsim_apply_dye_mixing(
    float *imageData,
    int width,
    int height,
    const filmsim_couplers_t *couplers
);

/* Film profile loading from JSON string */
bool filmsim_load_profile_from_json(
    filmsim_film_profile_t *profile,
    const char *json_string
);

/* Get embedded film profile by name */
const char* filmsim_get_embedded_profile(const char *film_name);

/* List all embedded film profiles */
int filmsim_list_embedded_profiles(char ***names);

#ifdef __cplusplus
}
#endif
