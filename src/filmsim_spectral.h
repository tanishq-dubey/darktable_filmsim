/*
 * Film Simulation Plugin for darktable
 * Spectral processing module
 *
 * RGB-to-spectral conversion using Jakob & Hanika 2019 method
 * Spectral exposure and dye formation models
 */

#pragma once

#include "filmsim_common.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Number of polynomial coefficients for RGB-to-spectral */
#define FILMSIM_SPECTRAL_N_COEFFS 3

/* Spectral range and resolution */
#define FILMSIM_SPECTRAL_START 400.0f  /* nm */
#define FILMSIM_SPECTRAL_END   700.0f  /* nm */
#define FILMSIM_SPECTRAL_STEP  10.0f   /* nm */
#define FILMSIM_SPECTRAL_BANDS 31      /* (700-400)/10 + 1 */

/* LUT dimensions for spectral processing */
#define FILMSIM_SPECTRAL_LUT_SIZE 32   /* 32^3 = 32768 entries */

/* RGB-to-spectral model data structure */
typedef struct {
    uint32_t res;      /* Resolution of the model (typically 64) */
    float *scale;      /* Scale values for z-axis lookup */
    float *data;       /* Coefficient data */
    bool loaded;       /* Whether model is successfully loaded */
} filmsim_rgb2spec_t;

/* Spectral LUT for fast processing */
typedef struct {
    /* Exposure LUT: RGB -> (E_Y, E_M, E_C) layer exposures */
    float *exposure_lut;    /* [LUT_SIZE^3 * 3] */
    
    /* Dye LUT: (D_Y, D_M, D_C) -> scanned RGB */
    float *dye_lut;         /* [LUT_SIZE^3 * 3] */
    
    /* Film base transmittance (from spectral_density_min) */
    float base_T_r, base_T_g, base_T_b;
    
    /* Kappa functions (dye absorption spectra) */
    float kappa_y[FILMSIM_SPECTRAL_BANDS];
    float kappa_m[FILMSIM_SPECTRAL_BANDS];
    float kappa_c[FILMSIM_SPECTRAL_BANDS];
    
    /* Film base density spectrum */
    float base_density[FILMSIM_SPECTRAL_BANDS];
    
    /* Film sensitivity spectra */
    float sens_y[FILMSIM_SPECTRAL_BANDS];
    float sens_m[FILMSIM_SPECTRAL_BANDS];
    float sens_c[FILMSIM_SPECTRAL_BANDS];
    
    bool valid;
} filmsim_spectral_lut_t;

/* Initialize/free RGB-to-spectral model */
bool filmsim_rgb2spec_init(filmsim_rgb2spec_t *model);
void filmsim_rgb2spec_free(filmsim_rgb2spec_t *model);

/* Load RGB-to-spectral coefficients from file */
bool filmsim_rgb2spec_load(filmsim_rgb2spec_t *model, const char *filename);

/* Fetch spectral coefficients for an RGB color */
void filmsim_rgb2spec_fetch(const filmsim_rgb2spec_t *model, 
                            const float rgb[3], 
                            float coeffs[FILMSIM_SPECTRAL_N_COEFFS]);

/* Evaluate spectrum at a wavelength given coefficients */
float filmsim_rgb2spec_eval(const float coeffs[FILMSIM_SPECTRAL_N_COEFFS], float lambda);

/* Build spectral LUTs for a film profile */
bool filmsim_spectral_build_luts(filmsim_spectral_lut_t *lut,
                                  const filmsim_film_profile_t *film,
                                  const filmsim_rgb2spec_t *rgb2spec);

/* Free spectral LUTs */
void filmsim_spectral_free_luts(filmsim_spectral_lut_t *lut);

/* Apply spectral exposure conversion (RGB -> layer exposures) */
void filmsim_spectral_apply_exposure(const filmsim_spectral_lut_t *lut,
                                      float *image_data,
                                      int width, int height);

/* Apply spectral dye model (dye densities -> scanned RGB) */
void filmsim_spectral_apply_dye(const filmsim_spectral_lut_t *lut,
                                 float *image_data,
                                 int width, int height);

/* Check if a film profile has spectral data */
bool filmsim_has_spectral_data(const filmsim_film_profile_t *film);

/* Interpolate spectral data at a wavelength */
float filmsim_interpolate_spectral(const filmsim_spectral_density_t *data, 
                                    int count, float wavelength);

#ifdef __cplusplus
}
#endif
