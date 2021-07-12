/**
 *  @file test_compress_ts.c
 *  @author Sheng Di
 *  @date May, 2018
 *  @brief This is an example of using compression interface
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <libgen.h>
#include <unistd.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include "tng_compress.h"

struct timeval startTime;
struct timeval endTime;  /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;
double elapsed = 0;


void cost_start() {
    totalCost = 0;
    gettimeofday(&costStart, NULL);
}

void cost_end() {
    struct timeval costEnd;
    gettimeofday(&costEnd, NULL);
    elapsed = ((costEnd.tv_sec * 1000000 + costEnd.tv_usec) - (costStart.tv_sec * 1000000 + costStart.tv_usec)) /
              1000000.0;
    totalCost += elapsed;
}

void verify(float *ori_data, float *data, size_t num_elements, double *psnr, double *nrmse, double *diffMax) {
    size_t i = 0;
    double Max = ori_data[0];
    double Min = ori_data[0];
    *diffMax = fabs(data[0] - ori_data[0]);
    double diff_sum = 0;
    double maxpw_relerr = fabs((data[0] - ori_data[0]) / ori_data[0]);
    double sum1 = 0, sum2 = 0;
    for (i = 0; i < num_elements; i++) {
        sum1 += ori_data[i];
        sum2 += data[i];
    }
    double mean1 = sum1 / num_elements;
    double mean2 = sum2 / num_elements;

    double sum3 = 0, sum4 = 0;
    double sum = 0, prodSum = 0, relerr = 0;

    double *diff = (double *) malloc(num_elements * sizeof(double));

    for (i = 0; i < num_elements; i++) {
        diff[i] = data[i] - ori_data[i];
        diff_sum += data[i] - ori_data[i];
        if (Max < ori_data[i]) Max = ori_data[i];
        if (Min > ori_data[i]) Min = ori_data[i];
        double err = fabs(data[i] - ori_data[i]);
        if (ori_data[i] != 0) {
            relerr = err / fabs(ori_data[i]);
            if (maxpw_relerr < relerr)
                maxpw_relerr = relerr;
        }

        if (*diffMax < err)
            *diffMax = err;
        prodSum += (ori_data[i] - mean1) * (data[i] - mean2);
        sum3 += (ori_data[i] - mean1) * (ori_data[i] - mean1);
        sum4 += (data[i] - mean2) * (data[i] - mean2);
        sum += err * err;
    }

    double mse = sum / num_elements;
    double range = Max - Min;
    *psnr = 20 * log10(range) - 10 * log10(mse);
    *nrmse = sqrt(mse) / range;

    free(diff);
}

void readfloat(char *filename, float *data, size_t nbEle) {
    FILE *pFile = fopen(filename, "rb");
    fread(data, 4, nbEle, pFile);
    fclose(pFile);
}

void writefloat(char *filename, float *data, size_t nbEle) {
    FILE *pFile = fopen(filename, "wb");
    fwrite(data, 4, nbEle, pFile);
    fclose(pFile);
}

int main(int argc, char *argv[]) {

    if (argc < 6) {
        printf("Test case: tng [filex] [filey] [filez] [ntimesteps] [natoms] [reb] blocksize \n");
        printf("Example: tng exaalt-5423x3137/x.f32.dat exaalt-5423x3137/y.f32.dat exaalt-5423x3137/z.f32.dat 5423 3137 1e-3 10\n");
        exit(0);
    }
    int argp = 1;
    char *filexname = argv[argp++];
    char *fileyname = argv[argp++];
    char *filezname = argv[argp++];
    size_t ntimesteps = atoi(argv[argp++]);
    size_t natoms = atoi(argv[argp++]);
    double reb = atof(argv[argp++]);
    int blocksize = atoi(argv[argp++]);

    size_t nbEle = ntimesteps * natoms;

    float *data_x = (float *) malloc(sizeof(float) * nbEle);
    float *data_y = (float *) malloc(sizeof(float) * nbEle);
    float *data_z = (float *) malloc(sizeof(float) * nbEle);
    float *dec_x = (float *) malloc(sizeof(float) * nbEle);
    float *dec_y = (float *) malloc(sizeof(float) * nbEle);
    float *dec_z = (float *) malloc(sizeof(float) * nbEle);
    float *data = (float *) malloc(sizeof(float) * blocksize * natoms * 3);
    float *dec_data = (float *) malloc(sizeof(float) * blocksize * natoms * 3);
    readfloat(filexname, data_x, nbEle);
    readfloat(fileyname, data_y, nbEle);
    readfloat(filezname, data_z, nbEle);

    int algo[4];
    algo[0] = -1;
    algo[1] = -1;
    algo[2] = -1;
    algo[3] = -1;
    size_t totalOutSize = 0;
    double totalCmprTime = 0, totalDecmprTime = 0;

    for (size_t ts = 0; ts < ntimesteps; ts += blocksize) {
        size_t block_timesteps = ts + blocksize > ntimesteps ? ntimesteps - ts : blocksize;
        size_t block_num = block_timesteps * natoms;
        int outSize;

        for (size_t i = 0; i < block_num; i++) {
            data[i * 3] = data_x[ts * natoms + i];
            data[i * 3 + 1] = data_y[ts * natoms + i];
            data[i * 3 + 2] = data_z[ts * natoms + i];
        }

        float min = data[0], max = data[0];
        for (size_t i = 0; i < block_num * 3; i++) {
            if (data[i] > max) {
                max = data[i];
            }
            if (data[i] < min) {
                min = data[i];
            }
        }
        float eb = reb * (max - min);
        if (ts == 0) {
            cost_start();
            int tuning_timesteps = block_timesteps > 5 ? 5 : block_timesteps;
            tng_compress_pos_float(data, natoms, tuning_timesteps, eb, 0, algo, &outSize);
            cost_end();
            totalCmprTime += elapsed;
            printf("Algo = %d %d %d %d\n", algo[0], algo[1], algo[2], algo[3]);
        }

        cost_start();
        char *compressed = tng_compress_pos_float(data, natoms, block_timesteps, eb, 0, algo, &outSize);
        cost_end();
        totalOutSize += outSize;
        totalCmprTime += elapsed;
//        printf("compressed_size=%d, compress_time=%.3f\n", outSize, elapsed);

        cost_start();
        tng_compress_uncompress_float(compressed, dec_data);
        cost_end();
        free(compressed);
        totalDecmprTime += elapsed;
//        printf("decompress_time=%.3f\n", elapsed);

        for (size_t i = 0; i < block_num; i++) {
            dec_x[ts * natoms + i] = dec_data[i * 3];
            dec_y[ts * natoms + i] = dec_data[i * 3 + 1];
            dec_z[ts * natoms + i] = dec_data[i * 3 + 2];
        }

    }

    double psnr, nrmse, max_diff;
    verify(data_x, dec_x, nbEle, &psnr, &nrmse, &max_diff);
    printf("method=tng, file=%s, block=%d, compression_ratio=%.3f, reb=%.1e, eb=%.6f, psnr=%.3f, nsmse=%e, compress_time=%.3f, decompress_time=%.3f\n",
           filexname, blocksize, nbEle * sizeof(float) * 3.0 / totalOutSize,
           reb, max_diff, psnr, nrmse, totalCmprTime / 3, totalDecmprTime / 3);

    verify(data_y, dec_y, nbEle, &psnr, &nrmse, &max_diff);
    printf("method=tng, file=%s, block=%d, compression_ratio=%.3f, reb=%.1e, eb=%.6f, psnr=%.3f, nsmse=%e, compress_time=%.3f, decompress_time=%.3f\n",
           fileyname, blocksize, nbEle * sizeof(float) * 3.0 / totalOutSize,
           reb, max_diff, psnr, nrmse, totalCmprTime / 3, totalDecmprTime / 3);

    verify(data_z, dec_z, nbEle, &psnr, &nrmse, &max_diff);
    printf("method=tng, file=%s, block=%d, compression_ratio=%.3f, reb=%.1e, eb=%.6f, psnr=%.3f, nsmse=%e, compress_time=%.3f, decompress_time=%.3f\n",
           filezname, blocksize, nbEle * sizeof(float) * 3.0 / totalOutSize,
           reb, max_diff, psnr, nrmse, totalCmprTime / 3, totalDecmprTime / 3);

    char outputname[600];
    sprintf(outputname, "%s.b%d.%.1e.tng.out",  basename(strdup(filexname)), blocksize, reb);
    writefloat(outputname, dec_x, nbEle);
    sprintf(outputname, "%s.b%d.%.1e.tng.out",  basename(strdup(fileyname)), blocksize, reb);
    writefloat(outputname, dec_y, nbEle);
    sprintf(outputname, "%s.b%d.%.1e.tng.out",  basename(strdup(filezname)), blocksize, reb);
    writefloat(outputname, dec_z, nbEle);

    free(data);
    free(dec_data);
    free(data_x);
    free(data_y);
    free(data_z);
    free(dec_x);
    free(dec_y);
    free(dec_z);

    return 0;
}
