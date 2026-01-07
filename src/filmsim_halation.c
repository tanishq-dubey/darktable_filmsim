/*
 * Film Simulation Plugin for darktable
 * Halation effect (light scatter from film base)
 *
 * Ported from src/Halation.cpp (C++)
 */

#include "filmsim_common.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>
#include <omp.h>

/* Next power of 2 for FFT-friendly sizes */
static int nextPow2(int n) {
    int p = 1;
    while (p < n) p *= 2;
    return p;
}

/* Downsample a single channel from interleaved RGB using box filter */
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

/* Generate exponential PSF kernel centered at origin for FFT convolution */
static void generateExponentialPSF(
    float *kernel, int fftW, int fftH, float sigma
) {
    memset(kernel, 0, fftW * fftH * sizeof(float));

    int extent = (int)ceilf(sigma * 5.0f);
    extent = (extent < fftW / 2) ? extent : fftW / 2;
    extent = (extent < fftH / 2) ? extent : fftH / 2;

    double sum = 0.0;

    for (int dy = -extent; dy <= extent; ++dy) {
        for (int dx = -extent; dx <= extent; ++dx) {
            float r = sqrtf((float)(dx * dx + dy * dy));
            float val = expf(-r / sigma);

            int y = (dy + fftH) % fftH;
            int x = (dx + fftW) % fftW;

            kernel[y * fftW + x] = val;
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

/* Perform FFT convolution with pre-computed kernel */
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

/* Bilinear sample from downsampled buffer */
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

void filmsim_apply_halation(
    float *imageData,
    int width,
    int height,
    const filmsim_halation_t *halation,
    float pixelSizeUm
) {
    float maxStrength = halation->strength.r;
    if (halation->strength.g > maxStrength) maxStrength = halation->strength.g;
    if (halation->strength.b > maxStrength) maxStrength = halation->strength.b;
    
    if (maxStrength < 0.001f) {
        return;
    }

    const int DS = 4;
    int dsW = (width + DS - 1) / DS;
    int dsH = (height + DS - 1) / DS;

    int fftW = nextPow2(dsW);
    int fftH = nextPow2(dsH);
    int fftN = fftW * fftH;

    float sigmaR = (float)(halation->size_um.r / pixelSizeUm / DS);
    float sigmaG = (float)(halation->size_um.g / pixelSizeUm / DS);
    float sigmaB = (float)(halation->size_um.b / pixelSizeUm / DS);

    float *chanR = (float*)calloc(fftN, sizeof(float));
    float *chanG = (float*)calloc(fftN, sizeof(float));
    float *chanB = (float*)calloc(fftN, sizeof(float));

    if (!chanR || !chanG || !chanB) {
        if (chanR) free(chanR);
        if (chanG) free(chanG);
        if (chanB) free(chanB);
        return;
    }

    downsampleChannel(imageData, width, height, chanR, dsW, dsH, fftW, fftH, 0, DS);
    downsampleChannel(imageData, width, height, chanG, dsW, dsH, fftW, fftH, 1, DS);
    downsampleChannel(imageData, width, height, chanB, dsW, dsH, fftW, fftH, 2, DS);

    int complexN = fftH * (fftW / 2 + 1);
    fftwf_complex *kerFreqR = fftwf_alloc_complex(complexN);
    fftwf_complex *kerFreqG = fftwf_alloc_complex(complexN);
    fftwf_complex *kerFreqB = fftwf_alloc_complex(complexN);

    float *psfR = (float*)calloc(fftN, sizeof(float));
    float *psfG = (float*)calloc(fftN, sizeof(float));
    float *psfB = (float*)calloc(fftN, sizeof(float));

    generateExponentialPSF(psfR, fftW, fftH, sigmaR);
    generateExponentialPSF(psfG, fftW, fftH, sigmaG);
    generateExponentialPSF(psfB, fftW, fftH, sigmaB);

    fftwf_plan fwdKerR = fftwf_plan_dft_r2c_2d(fftH, fftW, psfR, kerFreqR, FFTW_ESTIMATE);
    fftwf_plan fwdKerG = fftwf_plan_dft_r2c_2d(fftH, fftW, psfG, kerFreqG, FFTW_ESTIMATE);
    fftwf_plan fwdKerB = fftwf_plan_dft_r2c_2d(fftH, fftW, psfB, kerFreqB, FFTW_ESTIMATE);
    
    fftwf_execute(fwdKerR);
    fftwf_execute(fwdKerG);
    fftwf_execute(fwdKerB);
    
    fftwf_destroy_plan(fwdKerR);
    fftwf_destroy_plan(fwdKerG);
    fftwf_destroy_plan(fwdKerB);

    free(psfR);
    free(psfG);
    free(psfB);

    convolveFFTWithPrecomputedKernel(chanR, kerFreqR, fftW, fftH);
    convolveFFTWithPrecomputedKernel(chanG, kerFreqG, fftW, fftH);
    convolveFFTWithPrecomputedKernel(chanB, kerFreqB, fftW, fftH);

    fftwf_free(kerFreqR);
    fftwf_free(kerFreqG);
    fftwf_free(kerFreqB);

    float strengthR = (float)halation->strength.r;
    float strengthG = (float)halation->strength.g;
    float strengthB = (float)halation->strength.b;

    #pragma omp parallel for
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            int idx = (y * width + x) * 3;

            float dsX = (float)x / DS;
            float dsY = (float)y / DS;

            float haloR = bilinearSample(chanR, dsW, dsH, fftW, dsX, dsY);
            float haloG = bilinearSample(chanG, dsW, dsH, fftW, dsX, dsY);
            float haloB = bilinearSample(chanB, dsW, dsH, fftW, dsX, dsY);

            imageData[idx + 0] += strengthR * haloR;
            imageData[idx + 1] += strengthG * haloG;
            imageData[idx + 2] += strengthB * haloB;
        }
    }

    free(chanR);
    free(chanG);
    free(chanB);
}
