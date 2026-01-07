/*
 * Film Simulation Plugin for darktable
 * Spectral processing implementation
 *
 * RGB-to-spectral conversion using Jakob & Hanika 2019 method
 * Spectral exposure and dye formation models
 */

#include "filmsim_spectral.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <glib.h>
#include <glib/gstdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Helper functions */
static inline int filmsim_min_int(int a, int b) { return a < b ? a : b; }
static inline float filmsim_max_f(float a, float b) { return a > b ? a : b; }
static inline float filmsim_clamp_f(float x, float lo, float hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

/* FMA helper */
static inline float filmsim_fma(float a, float b, float c) {
    return a * b + c;
}

/*
 * Initialize RGB-to-spectral model
 */
bool filmsim_rgb2spec_init(filmsim_rgb2spec_t *model) {
    if (!model) return false;
    memset(model, 0, sizeof(filmsim_rgb2spec_t));
    return true;
}

/*
 * Free RGB-to-spectral model
 */
void filmsim_rgb2spec_free(filmsim_rgb2spec_t *model) {
    if (!model) return;
    if (model->scale) free(model->scale);
    if (model->data) free(model->data);
    model->scale = NULL;
    model->data = NULL;
    model->loaded = false;
}

/*
 * Load RGB-to-spectral coefficients from file
 * File format: "SPEC" magic, uint32 resolution, float[] scale, float[] data
 */
bool filmsim_rgb2spec_load(filmsim_rgb2spec_t *model, const char *filename) {
    if (!model || !filename) return false;
    
    filmsim_rgb2spec_free(model);
    
    FILE *f = g_fopen(filename, "rb");
    if (!f) {
        fprintf(stderr, "[filmsim] Cannot open spectral coefficients: %s\n", filename);
        return false;
    }
    
    char header[4];
    if (fread(header, 4, 1, f) != 1 || memcmp(header, "SPEC", 4) != 0) {
        fprintf(stderr, "[filmsim] Invalid spectral file header: %s\n", filename);
        fclose(f);
        return false;
    }
    
    if (fread(&model->res, sizeof(uint32_t), 1, f) != 1) {
        fprintf(stderr, "[filmsim] Cannot read resolution from: %s\n", filename);
        fclose(f);
        return false;
    }
    
    size_t size_scale = sizeof(float) * model->res;
    size_t size_data = sizeof(float) * model->res * model->res * 
                       model->res * 3 * FILMSIM_SPECTRAL_N_COEFFS;
    
    model->scale = (float *)malloc(size_scale);
    model->data = (float *)malloc(size_data);
    
    if (!model->scale || !model->data) {
        fprintf(stderr, "[filmsim] Out of memory loading spectral model\n");
        filmsim_rgb2spec_free(model);
        fclose(f);
        return false;
    }
    
    if (fread(model->scale, size_scale, 1, f) != 1 ||
        fread(model->data, size_data, 1, f) != 1) {
        fprintf(stderr, "[filmsim] Error reading spectral data from: %s\n", filename);
        filmsim_rgb2spec_free(model);
        fclose(f);
        return false;
    }
    
    fclose(f);
    model->loaded = true;
    fprintf(stderr, "[filmsim] Loaded spectral model (res=%u) from: %s\n", model->res, filename);
    return true;
}

/*
 * Find interval in sorted array (binary search)
 */
static int find_interval(const float *values, int size, float x) {
    int left = 0;
    int last_interval = size - 2;
    int len = last_interval;
    
    while (len > 0) {
        int half = len >> 1;
        int middle = left + half + 1;
        
        if (values[middle] < x) {
            left = middle;
            len -= half + 1;
        } else {
            len = half;
        }
    }
    
    return filmsim_min_int(left, last_interval);
}

/*
 * Fetch spectral coefficients for an RGB color
 * Uses trilinear interpolation in the coefficient table
 */
void filmsim_rgb2spec_fetch(const filmsim_rgb2spec_t *model, 
                            const float rgb[3], 
                            float coeffs[FILMSIM_SPECTRAL_N_COEFFS]) {
    if (!model || !model->loaded || !rgb || !coeffs) {
        /* Fallback: return flat spectrum */
        coeffs[0] = 0.0f;
        coeffs[1] = 0.0f;
        coeffs[2] = 0.5f;  /* Will evaluate to ~0.5 reflectance */
        return;
    }
    
    /* Clamp input RGB to [0, 1] */
    float r = filmsim_clamp_f(rgb[0], 0.0f, 1.0f);
    float g = filmsim_clamp_f(rgb[1], 0.0f, 1.0f);
    float b = filmsim_clamp_f(rgb[2], 0.0f, 1.0f);
    
    /* Determine largest RGB component */
    int i = 0;
    float rgb_clamped[3] = {r, g, b};
    for (int j = 1; j < 3; j++) {
        if (rgb_clamped[j] >= rgb_clamped[i]) {
            i = j;
        }
    }
    
    int res = (int)model->res;
    float z = rgb_clamped[i];
    
    /* Handle black/very dark colors */
    if (z < 1e-6f) {
        coeffs[0] = 0.0f;
        coeffs[1] = 0.0f;
        coeffs[2] = 0.0f;
        return;
    }
    
    float scale = (res - 1) / z;
    float x = rgb_clamped[(i + 1) % 3] * scale;
    float y = rgb_clamped[(i + 2) % 3] * scale;
    
    /* Bilinearly/trilinearly interpolated lookup */
    uint32_t xi = filmsim_min_int((uint32_t)x, res - 2);
    uint32_t yi = filmsim_min_int((uint32_t)y, res - 2);
    uint32_t zi = find_interval(model->scale, res, z);
    
    uint32_t offset = (((i * res + zi) * res + yi) * res + xi) * FILMSIM_SPECTRAL_N_COEFFS;
    uint32_t dx = FILMSIM_SPECTRAL_N_COEFFS;
    uint32_t dy = FILMSIM_SPECTRAL_N_COEFFS * res;
    uint32_t dz = FILMSIM_SPECTRAL_N_COEFFS * res * res;
    
    float x1 = x - xi, x0 = 1.0f - x1;
    float y1 = y - yi, y0 = 1.0f - y1;
    float z1 = (z - model->scale[zi]) / (model->scale[zi + 1] - model->scale[zi]);
    float z0 = 1.0f - z1;
    
    for (int j = 0; j < FILMSIM_SPECTRAL_N_COEFFS; j++) {
        coeffs[j] = ((model->data[offset               ] * x0 +
                      model->data[offset + dx          ] * x1) * y0 +
                     (model->data[offset + dy          ] * x0 +
                      model->data[offset + dy + dx     ] * x1) * y1) * z0 +
                    ((model->data[offset + dz          ] * x0 +
                      model->data[offset + dz + dx     ] * x1) * y0 +
                     (model->data[offset + dz + dy     ] * x0 +
                      model->data[offset + dz + dy + dx] * x1) * y1) * z1;
        offset++;
    }
}

/*
 * Evaluate spectrum at a wavelength given coefficients
 * lambda should be normalized to [0, 1] range where 0=400nm, 1=700nm
 */
float filmsim_rgb2spec_eval(const float coeffs[FILMSIM_SPECTRAL_N_COEFFS], float lambda) {
    /* Polynomial: x = c0*lambda^2 + c1*lambda + c2 */
    float x = filmsim_fma(filmsim_fma(coeffs[0], lambda, coeffs[1]), lambda, coeffs[2]);
    /* Sigmoid: 0.5 + 0.5 * x / sqrt(1 + x^2) */
    float y = 1.0f / sqrtf(filmsim_fma(x, x, 1.0f));
    return filmsim_fma(0.5f * x, y, 0.5f);
}

/*
 * Interpolate spectral data at a wavelength
 */
float filmsim_interpolate_spectral(const filmsim_spectral_density_t *data, 
                                    int count, float wavelength) {
    if (!data || count < 1) return 0.0f;
    
    /* Clamp to data range */
    if (wavelength <= data[0].wavelength) return (float)data[0].density;
    if (wavelength >= data[count-1].wavelength) return (float)data[count-1].density;
    
    /* Binary search for interval */
    int lo = 0, hi = count - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (data[mid].wavelength <= wavelength) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    
    /* Linear interpolation */
    float t = (wavelength - (float)data[lo].wavelength) / 
              ((float)data[hi].wavelength - (float)data[lo].wavelength);
    return (float)data[lo].density * (1.0f - t) + (float)data[hi].density * t;
}

/*
 * Interpolate spectral sensitivity at a wavelength (RGB version)
 */
static void interpolate_sensitivity_rgb(const filmsim_spectral_sensitivity_rgb_t *data,
                                         int count, float wavelength,
                                         float *out_y, float *out_m, float *out_c) {
    *out_y = *out_m = *out_c = 0.0f;
    if (!data || count < 1) return;
    
    if (wavelength <= data[0].wavelength) {
        *out_y = (float)data[0].y;
        *out_m = (float)data[0].m;
        *out_c = (float)data[0].c;
        return;
    }
    if (wavelength >= data[count-1].wavelength) {
        *out_y = (float)data[count-1].y;
        *out_m = (float)data[count-1].m;
        *out_c = (float)data[count-1].c;
        return;
    }
    
    /* Binary search */
    int lo = 0, hi = count - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (data[mid].wavelength <= wavelength) lo = mid;
        else hi = mid;
    }
    
    float t = (wavelength - (float)data[lo].wavelength) / 
              ((float)data[hi].wavelength - (float)data[lo].wavelength);
    *out_y = (float)data[lo].y * (1.0f - t) + (float)data[hi].y * t;
    *out_m = (float)data[lo].m * (1.0f - t) + (float)data[hi].m * t;
    *out_c = (float)data[lo].c * (1.0f - t) + (float)data[hi].c * t;
}

/*
 * Check if a film profile has spectral data
 */
bool filmsim_has_spectral_data(const filmsim_film_profile_t *film) {
    if (!film) return false;
    
    /* Need sensitivity data */
    if (film->num_spectral_sensitivity < 2) return false;
    
    /* Need either min/mid data (negatives) or y/m/c data (reversals) */
    bool has_min_mid = (film->num_spectral_density_min >= 2 && 
                        film->num_spectral_density_mid >= 2);
    bool has_ymc = (film->num_spectral_density_y >= 2 &&
                    film->num_spectral_density_m >= 2 &&
                    film->num_spectral_density_c >= 2);
    
    return has_min_mid || has_ymc;
}

/*
 * Build spectral LUTs for a film profile
 */
bool filmsim_spectral_build_luts(filmsim_spectral_lut_t *lut,
                                  const filmsim_film_profile_t *film,
                                  const filmsim_rgb2spec_t *rgb2spec) {
    if (!lut || !film) return false;
    
    memset(lut, 0, sizeof(filmsim_spectral_lut_t));
    
    /* Check if we have spectral data */
    if (!filmsim_has_spectral_data(film)) {
        fprintf(stderr, "[filmsim] Film profile lacks spectral data\n");
        return false;
    }
    
    /* Allocate LUTs */
    size_t lut_entries = FILMSIM_SPECTRAL_LUT_SIZE * FILMSIM_SPECTRAL_LUT_SIZE * 
                         FILMSIM_SPECTRAL_LUT_SIZE;
    lut->exposure_lut = (float *)calloc(lut_entries * 3, sizeof(float));
    lut->dye_lut = (float *)calloc(lut_entries * 3, sizeof(float));
    
    if (!lut->exposure_lut || !lut->dye_lut) {
        filmsim_spectral_free_luts(lut);
        return false;
    }
    
    /* Build wavelength samples */
    float wavelengths[FILMSIM_SPECTRAL_BANDS];
    for (int i = 0; i < FILMSIM_SPECTRAL_BANDS; i++) {
        wavelengths[i] = FILMSIM_SPECTRAL_START + i * FILMSIM_SPECTRAL_STEP;
    }
    
    /* Interpolate film sensitivity and dye spectra */
    for (int i = 0; i < FILMSIM_SPECTRAL_BANDS; i++) {
        float wl = wavelengths[i];
        
        /* Sensitivity (Y, M, C layers) */
        interpolate_sensitivity_rgb(film->spectral_sensitivity.rgb,
                                     film->num_spectral_sensitivity, wl,
                                     &lut->sens_y[i], &lut->sens_m[i], &lut->sens_c[i]);
        
        /* Base density (from spectral_density_min) */
        if (film->num_spectral_density_min > 0) {
            lut->base_density[i] = filmsim_interpolate_spectral(
                film->spectral_density_min, film->num_spectral_density_min, wl);
        }
        
        /* Kappa functions (dye absorption spectra) */
        if (film->num_spectral_density_y > 0) {
            /* Reversal film: direct Y/M/C dye data */
            lut->kappa_y[i] = filmsim_interpolate_spectral(
                film->spectral_density_y, film->num_spectral_density_y, wl);
            lut->kappa_m[i] = filmsim_interpolate_spectral(
                film->spectral_density_m, film->num_spectral_density_m, wl);
            lut->kappa_c[i] = filmsim_interpolate_spectral(
                film->spectral_density_c, film->num_spectral_density_c, wl);
        } else if (film->num_spectral_density_min > 0 && film->num_spectral_density_mid > 0) {
            /* Negative film: derive from min/mid difference */
            float d_min = lut->base_density[i];
            float d_mid = filmsim_interpolate_spectral(
                film->spectral_density_mid, film->num_spectral_density_mid, wl);
            float d_dye = filmsim_max_f(0.0f, d_mid - d_min);
            
            /* Simple wavelength-based partitioning into Y/M/C */
            /* Y peaks at ~450nm, M at ~550nm, C at ~650nm */
            float sigma = 30.0f;
            float g_y = expf(-0.5f * powf((wl - 450.0f) / sigma, 2.0f));
            float g_m = expf(-0.5f * powf((wl - 550.0f) / sigma, 2.0f));
            float g_c = expf(-0.5f * powf((wl - 650.0f) / sigma, 2.0f));
            float sum = g_y + g_m + g_c + 1e-6f;
            
            lut->kappa_y[i] = d_dye * g_y / sum;
            lut->kappa_m[i] = d_dye * g_m / sum;
            lut->kappa_c[i] = d_dye * g_c / sum;
        }
    }
    
    /* Normalize kappa functions to peak = 1.0 */
    float max_y = 0.0f, max_m = 0.0f, max_c = 0.0f;
    for (int i = 0; i < FILMSIM_SPECTRAL_BANDS; i++) {
        if (lut->kappa_y[i] > max_y) max_y = lut->kappa_y[i];
        if (lut->kappa_m[i] > max_m) max_m = lut->kappa_m[i];
        if (lut->kappa_c[i] > max_c) max_c = lut->kappa_c[i];
    }
    if (max_y > 1e-6f) for (int i = 0; i < FILMSIM_SPECTRAL_BANDS; i++) lut->kappa_y[i] /= max_y;
    if (max_m > 1e-6f) for (int i = 0; i < FILMSIM_SPECTRAL_BANDS; i++) lut->kappa_m[i] /= max_m;
    if (max_c > 1e-6f) for (int i = 0; i < FILMSIM_SPECTRAL_BANDS; i++) lut->kappa_c[i] /= max_c;
    
    /* Compute base transmittance (average over spectral regions) */
    float base_r = 0.0f, base_g = 0.0f, base_b = 0.0f;
    int count_r = 0, count_g = 0, count_b = 0;
    for (int i = 0; i < FILMSIM_SPECTRAL_BANDS; i++) {
        float wl = wavelengths[i];
        float d = lut->base_density[i];
        if (wl >= 600) { base_r += d; count_r++; }
        else if (wl >= 500) { base_g += d; count_g++; }
        else { base_b += d; count_b++; }
    }
    if (count_r > 0) base_r /= count_r;
    if (count_g > 0) base_g /= count_g;
    if (count_b > 0) base_b /= count_b;
    
    lut->base_T_r = powf(10.0f, -base_r);
    lut->base_T_g = powf(10.0f, -base_g);
    lut->base_T_b = powf(10.0f, -base_b);
    
    /* Build exposure LUT: RGB -> layer exposures */
    /* This integrates input spectrum against film sensitivity curves */
    const int LUT_SIZE = FILMSIM_SPECTRAL_LUT_SIZE;
    
    #ifdef _OPENMP
    #pragma omp parallel for collapse(3)
    #endif
    for (int ir = 0; ir < LUT_SIZE; ir++) {
        for (int ig = 0; ig < LUT_SIZE; ig++) {
            for (int ib = 0; ib < LUT_SIZE; ib++) {
                float r = (float)ir / (LUT_SIZE - 1);
                float g = (float)ig / (LUT_SIZE - 1);
                float b = (float)ib / (LUT_SIZE - 1);
                
                /* Convert RGB to spectral coefficients */
                float rgb[3] = {r, g, b};
                float coeffs[FILMSIM_SPECTRAL_N_COEFFS];
                
                if (rgb2spec && rgb2spec->loaded) {
                    filmsim_rgb2spec_fetch(rgb2spec, rgb, coeffs);
                } else {
                    /* Fallback: approximate spectrum from RGB */
                    coeffs[0] = 0.0f;
                    coeffs[1] = 0.0f;
                    coeffs[2] = 0.0f;
                }
                
                /* Integrate spectrum against sensitivity curves */
                float E_y = 0.0f, E_m = 0.0f, E_c = 0.0f;
                for (int wi = 0; wi < FILMSIM_SPECTRAL_BANDS; wi++) {
                    float lambda_norm = (float)wi / (FILMSIM_SPECTRAL_BANDS - 1);
                    float spectrum;
                    
                    if (rgb2spec && rgb2spec->loaded) {
                        spectrum = filmsim_rgb2spec_eval(coeffs, lambda_norm);
                    } else {
                        /* Simple RGB to spectrum approximation */
                        float wl = wavelengths[wi];
                        if (wl < 500) spectrum = b;
                        else if (wl < 600) spectrum = g;
                        else spectrum = r;
                    }
                    
                    E_y += spectrum * lut->sens_y[wi];
                    E_m += spectrum * lut->sens_m[wi];
                    E_c += spectrum * lut->sens_c[wi];
                }
                
                /* Normalize by number of bands */
                E_y /= FILMSIM_SPECTRAL_BANDS;
                E_m /= FILMSIM_SPECTRAL_BANDS;
                E_c /= FILMSIM_SPECTRAL_BANDS;
                
                /* Store in LUT */
                size_t idx = ((ir * LUT_SIZE + ig) * LUT_SIZE + ib) * 3;
                lut->exposure_lut[idx + 0] = E_y;
                lut->exposure_lut[idx + 1] = E_m;
                lut->exposure_lut[idx + 2] = E_c;
            }
        }
    }
    
    /* Build dye LUT: dye densities -> scanned RGB */
    /* Uses 5000K illuminant and standard scanner sensitivities */
    float illuminant[FILMSIM_SPECTRAL_BANDS];
    float scanner_r[FILMSIM_SPECTRAL_BANDS];
    float scanner_g[FILMSIM_SPECTRAL_BANDS];
    float scanner_b[FILMSIM_SPECTRAL_BANDS];
    
    /* Planckian illuminant at 5000K */
    const float T = 5000.0f;
    for (int i = 0; i < FILMSIM_SPECTRAL_BANDS; i++) {
        float wl = wavelengths[i] * 1e-9f;  /* Convert nm to m */
        float h = 6.626e-34f;
        float c = 3e8f;
        float k = 1.381e-23f;
        illuminant[i] = (2.0f * h * c * c) / (powf(wl, 5.0f) * (expf(h * c / (wl * k * T)) - 1.0f));
        
        /* Gaussian scanner sensitivities */
        float wl_nm = wavelengths[i];
        scanner_r[i] = expf(-0.5f * powf((wl_nm - 670.0f) / 53.0f, 2.0f));
        scanner_g[i] = expf(-0.5f * powf((wl_nm - 532.0f) / 53.0f, 2.0f));
        scanner_b[i] = expf(-0.5f * powf((wl_nm - 473.0f) / 53.0f, 2.0f));
    }
    
    /* Normalize illuminant */
    float ill_sum = 0.0f;
    for (int i = 0; i < FILMSIM_SPECTRAL_BANDS; i++) ill_sum += illuminant[i];
    if (ill_sum > 0) for (int i = 0; i < FILMSIM_SPECTRAL_BANDS; i++) illuminant[i] /= ill_sum;
    
    #ifdef _OPENMP
    #pragma omp parallel for collapse(3)
    #endif
    for (int iy = 0; iy < LUT_SIZE; iy++) {
        for (int im = 0; im < LUT_SIZE; im++) {
            for (int ic = 0; ic < LUT_SIZE; ic++) {
                /* Dye densities in range [0, 3] */
                float D_y = (float)iy / (LUT_SIZE - 1) * 3.0f;
                float D_m = (float)im / (LUT_SIZE - 1) * 3.0f;
                float D_c = (float)ic / (LUT_SIZE - 1) * 3.0f;
                
                /* Integrate spectral transmission against scanner */
                float scan_r = 0.0f, scan_g = 0.0f, scan_b = 0.0f;
                for (int wi = 0; wi < FILMSIM_SPECTRAL_BANDS; wi++) {
                    /* Total density = base + dye contributions */
                    float D_total = lut->base_density[wi] +
                                    D_y * lut->kappa_y[wi] +
                                    D_m * lut->kappa_m[wi] +
                                    D_c * lut->kappa_c[wi];
                    
                    /* Transmittance */
                    float T_wl = powf(10.0f, -D_total);
                    
                    /* Scanner response */
                    float signal = illuminant[wi] * T_wl;
                    scan_r += signal * scanner_r[wi];
                    scan_g += signal * scanner_g[wi];
                    scan_b += signal * scanner_b[wi];
                }
                
                /* Normalize */
                scan_r /= FILMSIM_SPECTRAL_BANDS;
                scan_g /= FILMSIM_SPECTRAL_BANDS;
                scan_b /= FILMSIM_SPECTRAL_BANDS;
                
                /* Store in LUT */
                size_t idx = ((iy * LUT_SIZE + im) * LUT_SIZE + ic) * 3;
                lut->dye_lut[idx + 0] = scan_r;
                lut->dye_lut[idx + 1] = scan_g;
                lut->dye_lut[idx + 2] = scan_b;
            }
        }
    }
    
    lut->valid = true;
    fprintf(stderr, "[filmsim] Built spectral LUTs (%dx%dx%d)\n", 
            LUT_SIZE, LUT_SIZE, LUT_SIZE);
    return true;
}

/*
 * Free spectral LUTs
 */
void filmsim_spectral_free_luts(filmsim_spectral_lut_t *lut) {
    if (!lut) return;
    if (lut->exposure_lut) free(lut->exposure_lut);
    if (lut->dye_lut) free(lut->dye_lut);
    memset(lut, 0, sizeof(filmsim_spectral_lut_t));
}

/*
 * Trilinear interpolation in 3D LUT
 */
static void lut_trilinear(const float *lut, int size, 
                          float r, float g, float b,
                          float scale,  /* Input range scale factor */
                          float *out_r, float *out_g, float *out_b) {
    /* Scale inputs to LUT range */
    float fr = filmsim_clamp_f(r / scale, 0.0f, 1.0f) * (size - 1);
    float fg = filmsim_clamp_f(g / scale, 0.0f, 1.0f) * (size - 1);
    float fb = filmsim_clamp_f(b / scale, 0.0f, 1.0f) * (size - 1);
    
    int ir0 = filmsim_min_int((int)fr, size - 2);
    int ig0 = filmsim_min_int((int)fg, size - 2);
    int ib0 = filmsim_min_int((int)fb, size - 2);
    
    float dr = fr - ir0;
    float dg = fg - ig0;
    float db = fb - ib0;
    
    /* 8 corner lookups */
    #define LUT_IDX(ir, ig, ib) (((ir) * size + (ig)) * size + (ib)) * 3
    
    size_t i000 = LUT_IDX(ir0, ig0, ib0);
    size_t i001 = LUT_IDX(ir0, ig0, ib0+1);
    size_t i010 = LUT_IDX(ir0, ig0+1, ib0);
    size_t i011 = LUT_IDX(ir0, ig0+1, ib0+1);
    size_t i100 = LUT_IDX(ir0+1, ig0, ib0);
    size_t i101 = LUT_IDX(ir0+1, ig0, ib0+1);
    size_t i110 = LUT_IDX(ir0+1, ig0+1, ib0);
    size_t i111 = LUT_IDX(ir0+1, ig0+1, ib0+1);
    
    #undef LUT_IDX
    
    /* Trilinear interpolation for each output channel */
    for (int c = 0; c < 3; c++) {
        float c00 = lut[i000+c] * (1-dr) + lut[i100+c] * dr;
        float c01 = lut[i001+c] * (1-dr) + lut[i101+c] * dr;
        float c10 = lut[i010+c] * (1-dr) + lut[i110+c] * dr;
        float c11 = lut[i011+c] * (1-dr) + lut[i111+c] * dr;
        
        float c0 = c00 * (1-dg) + c10 * dg;
        float c1 = c01 * (1-dg) + c11 * dg;
        
        float result = c0 * (1-db) + c1 * db;
        
        if (c == 0) *out_r = result;
        else if (c == 1) *out_g = result;
        else *out_b = result;
    }
}

/*
 * Apply spectral exposure conversion (RGB -> layer exposures)
 */
void filmsim_spectral_apply_exposure(const filmsim_spectral_lut_t *lut,
                                      float *image_data,
                                      int width, int height) {
    if (!lut || !lut->valid || !lut->exposure_lut || !image_data) return;
    
    const size_t npixels = (size_t)width * height;
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_t i = 0; i < npixels; i++) {
        float r = image_data[i * 3 + 0];
        float g = image_data[i * 3 + 1];
        float b = image_data[i * 3 + 2];
        
        float E_y, E_m, E_c;
        lut_trilinear(lut->exposure_lut, FILMSIM_SPECTRAL_LUT_SIZE,
                      r, g, b, 1.0f,
                      &E_y, &E_m, &E_c);
        
        /* Output: layer exposures (will be processed by H-D curves) */
        /* Map back to R/G/B channels for curve processing:
         * R (red-sensitive layer) -> Cyan dye -> E_c
         * G (green-sensitive layer) -> Magenta dye -> E_m
         * B (blue-sensitive layer) -> Yellow dye -> E_y
         */
        image_data[i * 3 + 0] = E_c;
        image_data[i * 3 + 1] = E_m;
        image_data[i * 3 + 2] = E_y;
    }
}

/*
 * Apply spectral dye model (dye densities -> scanned RGB)
 */
void filmsim_spectral_apply_dye(const filmsim_spectral_lut_t *lut,
                                 float *image_data,
                                 int width, int height) {
    if (!lut || !lut->valid || !lut->dye_lut || !image_data) return;
    
    const size_t npixels = (size_t)width * height;
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_t i = 0; i < npixels; i++) {
        /* Input: D_C, D_M, D_Y (from H-D curves + dye mixing) */
        float D_c = image_data[i * 3 + 0];
        float D_m = image_data[i * 3 + 1];
        float D_y = image_data[i * 3 + 2];
        
        float scan_r, scan_g, scan_b;
        lut_trilinear(lut->dye_lut, FILMSIM_SPECTRAL_LUT_SIZE,
                      D_y, D_m, D_c, 3.0f,  /* Density range 0-3 */
                      &scan_r, &scan_g, &scan_b);
        
        /* Normalize by film base transmittance */
        if (lut->base_T_r > 1e-6f) scan_r /= lut->base_T_r;
        if (lut->base_T_g > 1e-6f) scan_g /= lut->base_T_g;
        if (lut->base_T_b > 1e-6f) scan_b /= lut->base_T_b;
        
        image_data[i * 3 + 0] = scan_r;
        image_data[i * 3 + 1] = scan_g;
        image_data[i * 3 + 2] = scan_b;
    }
}
