/*
 * Film Simulation Plugin for darktable
 * Monotone cubic spline interpolation
 *
 * Ported from src/Spline.h (C++ template class)
 */

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
    if (n < 2) {
        spline->x = NULL;
        spline->y = NULL;
        spline->m = NULL;
        spline->n = 0;
        return;
    }

    /* Allocate memory */
    spline->n = n;
    spline->x = (double*)malloc(n * sizeof(double));
    spline->y = (double*)malloc(n * sizeof(double));
    spline->m = (double*)malloc(n * sizeof(double));

    if (!spline->x || !spline->y || !spline->m) {
        /* Memory allocation failed */
        free(spline->x);
        free(spline->y);
        free(spline->m);
        spline->x = spline->y = spline->m = NULL;
        spline->n = 0;
        return;
    }

    /* Copy data */
    memcpy(spline->x, x, n * sizeof(double));
    memcpy(spline->y, y, n * sizeof(double));

    /* Compute secants */
    double *d = (double*)malloc((n-1) * sizeof(double));
    if (!d) {
        free(spline->x);
        free(spline->y);
        free(spline->m);
        spline->x = spline->y = spline->m = NULL;
        spline->n = 0;
        return;
    }

    for (int i = 0; i < n-1; i++) {
        double h = spline->x[i+1] - spline->x[i];
        if (h <= 0) {
            /* Invalid or unsorted data */
            free(d);
            free(spline->x);
            free(spline->y);
            free(spline->m);
            spline->x = spline->y = spline->m = NULL;
            spline->n = 0;
            return;
        }
        d[i] = (spline->y[i+1] - spline->y[i]) / h;
    }

    /* Initialize tangents (Fritsch-Carlson algorithm) */
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

    /* Ensure monotonicity (Fritsch-Carlson) */
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

    /* Binary search for segment */
    int i = 0;
    int j = spline->n - 1;
    while (i < j) {
        int m = (i + j) / 2;
        if (x < spline->x[m]) {
            j = m;
        } else {
            i = m + 1;
        }
    }
    i--;

    /* Hermite interpolation */
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
    if (spline->x) free(spline->x);
    if (spline->y) free(spline->y);
    if (spline->m) free(spline->m);
    spline->x = spline->y = spline->m = NULL;
    spline->n = 0;
}
