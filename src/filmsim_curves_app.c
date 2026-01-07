/*
 * Film Simulation Plugin for darktable
 * Film profile application (H-D curves) and dye mixing
 */

#include "filmsim_common.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define LUT_SIZE 65536
#define LN10 2.302585092994046

void filmsim_apply_dye_mixing(
    float *imageData,
    int width,
    int height,
    const filmsim_couplers_t *couplers
) {
    const double *mix = (double*)couplers->dye_mix_matrix;

    #pragma omp parallel for
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            int idx = (y * width + x) * 3;

            float D_C = imageData[idx];
            float D_M = imageData[idx + 1];
            float D_Y = imageData[idx + 2];

            float C_out = (float)(mix[0] * D_C + mix[1] * D_M + mix[2] * D_Y);
            float M_out = (float)(mix[3] * D_C + mix[4] * D_M + mix[5] * D_Y);
            float Y_out = (float)(mix[6] * D_C + mix[7] * D_M + mix[8] * D_Y);

            imageData[idx]     = C_out;
            imageData[idx + 1] = M_out;
            imageData[idx + 2] = Y_out;
        }
    }
}

void filmsim_apply_film_profile(
    float *imageData,
    int width,
    int height,
    const filmsim_film_profile_t *film,
    float pushPullStops
) {
    filmsim_spline_t splineR, splineG, splineB;

    double xvals[FILMSIM_MAX_CURVE_POINTS];
    double yR[FILMSIM_MAX_CURVE_POINTS];
    double yG[FILMSIM_MAX_CURVE_POINTS];
    double yB[FILMSIM_MAX_CURVE_POINTS];

    int n = film->num_curves;
    if (n == 0 || n > FILMSIM_MAX_CURVE_POINTS) return;

    for (int i = 0; i < n; i++) {
        xvals[i] = film->curves.rgb[i].d;
        yR[i] = film->curves.rgb[i].r;
        yG[i] = film->curves.rgb[i].g;
        yB[i] = film->curves.rgb[i].b;
    }

    filmsim_spline_fit(&splineR, xvals, yR, n);
    filmsim_spline_fit(&splineG, xvals, yG, n);
    filmsim_spline_fit(&splineB, xvals, yB, n);

    float lutR[LUT_SIZE], lutG[LUT_SIZE], lutB[LUT_SIZE];
    float push_pull_log = pushPullStops * (float)log10(2.0);

    for (int i = 0; i < LUT_SIZE; i++) {
        float exposure = (float)i / (LUT_SIZE - 1.0f);

        if (exposure < 1e-8f) exposure = 1e-8f;
        float logE = log10f(exposure) + push_pull_log;

        lutR[i] = (float)filmsim_spline_interpolate(&splineR, logE);
        lutG[i] = (float)filmsim_spline_interpolate(&splineG, logE);
        lutB[i] = (float)filmsim_spline_interpolate(&splineB, logE);
    }

    #pragma omp parallel for
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            int idx = (y * width + x) * 3;

            float r = imageData[idx];
            float g = imageData[idx + 1];
            float b = imageData[idx + 2];

            int idxR = (int)(r * (LUT_SIZE - 1));
            int idxG = (int)(g * (LUT_SIZE - 1));
            int idxB = (int)(b * (LUT_SIZE - 1));

            if (idxR < 0) idxR = 0;
            if (idxR >= LUT_SIZE) idxR = LUT_SIZE - 1;
            if (idxG < 0) idxG = 0;
            if (idxG >= LUT_SIZE) idxG = LUT_SIZE - 1;
            if (idxB < 0) idxB = 0;
            if (idxB >= LUT_SIZE) idxB = LUT_SIZE - 1;

            imageData[idx]     = lutR[idxR];
            imageData[idx + 1] = lutG[idxG];
            imageData[idx + 2] = lutB[idxB];
        }
    }

    filmsim_spline_free(&splineR);
    filmsim_spline_free(&splineG);
    filmsim_spline_free(&splineB);
}
