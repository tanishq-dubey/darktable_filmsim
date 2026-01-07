/*
 * Film Simulation Plugin for darktable
 * DIR (Development Inhibitor Releasing) edge enhancement
 *
 * Ported from src/DIR.cpp (C++)
 */

#include "filmsim_common.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>
#include <omp.h>

static int nextPow2(int n) {
    int p = 1;
    while (p < n) p *= 2;
    return p;
}

static void downsampleChannel(
    const float *src, int srcW, int srcH,
    float *dst, int dstW, int dstH, int fftW, int fftH,
    int channel, int downsampleFactor
) {
    const int DS = downsampleFactor;

    #pragma omp parallel for
    for (int dy = 0; dy < dstH; ++dy) {
        for (int dx = 0; dx < dstW; ++dx) {
            float sum = 0.0f;
            int count = 0;

            int sy0 = dy * DS;
            int sx0 = dx * DS;
            int sy1 = (sy0 + DS < srcH) ? sy0 + DS : srcH;
            int sx1 = (sx0 + DS < srcW) ? sx0 + DS : srcW;

            for (int sy = sy0; sy < sy1; ++sy) {
                for (int sx = sx0; sx < sx1; ++sx) {
                    sum += src[(sy * srcW + sx) * 3 + channel];
                    ++count;
                }
            }

            dst[dy * fftW + dx] = (count > 0) ? sum / count : 0.0f;
        }
    }
}

static void generateGaussianPSF(
    float *kernel, int fftW, int fftH, float sigma
) {
    memset(kernel, 0, fftW * fftH * sizeof(float));

    int extent = (int)ceilf(sigma * 3.0f);
    extent = (extent < fftW / 2) ? extent : fftW / 2;
    extent = (extent < fftH / 2) ? extent : fftH / 2;

    double sum = 0.0;
    double twoSigmaSq = 2.0 * sigma * sigma;

    for (int dy = -extent; dy <= extent; ++dy) {
        for (int dx = -extent; dx <= extent; ++dx) {
            double rSq = (double)(dx * dx + dy * dy);
            double val = exp(-rSq / twoSigmaSq);

            int y = (dy + fftH) % fftH;
            int x = (dx + fftW) % fftW;

            kernel[y * fftW + x] = (float)val;
            sum += val;
        }
    }

    if (sum > 0) {
        float invSum = (float)(1.0 / sum);
        for (int i = 0; i < fftW * fftH; ++i) {
            kernel[i] *= invSum;
        }
    }
}

static void convolveFFTWithPrecomputedKernel(
    float *image,
    const fftwf_complex *kerFreq,
    int fftW, int fftH
) {
    int n = fftW * fftH;
    int complexN = fftH * (fftW / 2 + 1);

    fftwf_complex *imgFreq = fftwf_alloc_complex(complexN);

    fftwf_plan fwdImg = fftwf_plan_dft_r2c_2d(fftH, fftW, image, imgFreq, FFTW_ESTIMATE);
    fftwf_plan inv = fftwf_plan_dft_c2r_2d(fftH, fftW, imgFreq, image, FFTW_ESTIMATE);

    fftwf_execute(fwdImg);

    #pragma omp parallel for
    for (int i = 0; i < complexN; ++i) {
        float a = imgFreq[i][0];
        float b = imgFreq[i][1];
        float c = kerFreq[i][0];
        float d = kerFreq[i][1];

        imgFreq[i][0] = a * c - b * d;
        imgFreq[i][1] = a * d + b * c;
    }

    fftwf_execute(inv);

    float scale = 1.0f / n;
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        image[i] *= scale;
    }

    fftwf_destroy_plan(fwdImg);
    fftwf_destroy_plan(inv);
    fftwf_free(imgFreq);
}

static float bilinearSample(
    const float *buf, int bufW, int bufH, int fftW,
    float x, float y
) {
    if (x < 0.0f) x = 0.0f;
    if (x > bufW - 1.0f) x = (float)(bufW - 1);
    if (y < 0.0f) y = 0.0f;
    if (y > bufH - 1.0f) y = (float)(bufH - 1);

    int x0 = (int)x;
    int y0 = (int)y;
    int x1 = (x0 + 1 < bufW) ? x0 + 1 : bufW - 1;
    int y1 = (y0 + 1 < bufH) ? y0 + 1 : bufH - 1;

    float fx = x - (float)x0;
    float fy = y - (float)y0;

    float v00 = buf[y0 * fftW + x0];
    float v10 = buf[y0 * fftW + x1];
    float v01 = buf[y1 * fftW + x0];
    float v11 = buf[y1 * fftW + x1];

    return v00 * (1.0f - fx) * (1.0f - fy) +
           v10 * fx * (1.0f - fy) +
           v01 * (1.0f - fx) * fy +
           v11 * fx * fy;
}

void filmsim_apply_dir(
    float *imageData,
    int width,
    int height,
    const filmsim_couplers_t *couplers,
    float pixelSizeUm
) {
    if (couplers->dir_amount_rgb[0] < 0.001f &&
        couplers->dir_amount_rgb[1] < 0.001f &&
        couplers->dir_amount_rgb[2] < 0.001f) {
        return;
    }

    const int DS = 2;
    int dsW = (width + DS - 1) / DS;
    int dsH = (height + DS - 1) / DS;

    int fftW = nextPow2(dsW);
    int fftH = nextPow2(dsH);
    int fftN = fftW * fftH;

    float sigma = (float)(couplers->dir_diffusion_um / pixelSizeUm / DS);

    float *chanC = (float*)calloc(fftN, sizeof(float));
    float *chanM = (float*)calloc(fftN, sizeof(float));
    float *chanY = (float*)calloc(fftN, sizeof(float));

    float *origC = (float*)calloc(fftN, sizeof(float));
    float *origM = (float*)calloc(fftN, sizeof(float));
    float *origY = (float*)calloc(fftN, sizeof(float));

    if (!chanC || !chanM || !chanY || !origC || !origM || !origY) {
        if (chanC) free(chanC);
        if (chanM) free(chanM);
        if (chanY) free(chanY);
        if (origC) free(origC);
        if (origM) free(origM);
        if (origY) free(origY);
        return;
    }

    downsampleChannel(imageData, width, height, chanC, dsW, dsH, fftW, fftH, 0, DS);
    downsampleChannel(imageData, width, height, chanM, dsW, dsH, fftW, fftH, 1, DS);
    downsampleChannel(imageData, width, height, chanY, dsW, dsH, fftW, fftH, 2, DS);

    memcpy(origC, chanC, fftN * sizeof(float));
    memcpy(origM, chanM, fftN * sizeof(float));
    memcpy(origY, chanY, fftN * sizeof(float));

    int complexN = fftH * (fftW / 2 + 1);
    fftwf_complex *kerFreq = fftwf_alloc_complex(complexN);

    float *psf = (float*)calloc(fftN, sizeof(float));
    generateGaussianPSF(psf, fftW, fftH, sigma);

    fftwf_plan fwdKer = fftwf_plan_dft_r2c_2d(fftH, fftW, psf, kerFreq, FFTW_ESTIMATE);
    fftwf_execute(fwdKer);
    fftwf_destroy_plan(fwdKer);

    free(psf);

    convolveFFTWithPrecomputedKernel(chanC, kerFreq, fftW, fftH);
    convolveFFTWithPrecomputedKernel(chanM, kerFreq, fftW, fftH);
    convolveFFTWithPrecomputedKernel(chanY, kerFreq, fftW, fftH);

    fftwf_free(kerFreq);

    float amountC = (float)couplers->dir_amount_rgb[0];
    float amountM = (float)couplers->dir_amount_rgb[1];
    float amountY = (float)couplers->dir_amount_rgb[2];

    #pragma omp parallel for
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            int idx = (y * width + x) * 3;

            float dsX = (float)x / DS;
            float dsY = (float)y / DS;

            float origCVal = bilinearSample(origC, dsW, dsH, fftW, dsX, dsY);
            float origMVal = bilinearSample(origM, dsW, dsH, fftW, dsX, dsY);
            float origYVal = bilinearSample(origY, dsW, dsH, fftW, dsX, dsY);

            float blurCVal = bilinearSample(chanC, dsW, dsH, fftW, dsX, dsY);
            float blurMVal = bilinearSample(chanM, dsW, dsH, fftW, dsX, dsY);
            float blurYVal = bilinearSample(chanY, dsW, dsH, fftW, dsX, dsY);

            float D_C = imageData[idx];
            float D_M = imageData[idx + 1];
            float D_Y = imageData[idx + 2];

            float diffC = origCVal - blurCVal;
            float diffM = origMVal - blurMVal;
            float diffY = origYVal - blurYVal;

            imageData[idx]     = D_C + amountC * diffC;
            imageData[idx + 1] = D_M + amountM * diffM;
            imageData[idx + 2] = D_Y + amountY * diffY;

            if (imageData[idx] < 0.0f) imageData[idx] = 0.0f;
            if (imageData[idx + 1] < 0.0f) imageData[idx + 1] = 0.0f;
            if (imageData[idx + 2] < 0.0f) imageData[idx + 2] = 0.0f;
        }
    }

    free(chanC);
    free(chanM);
    free(chanY);
    free(origC);
    free(origM);
    free(origY);
}
