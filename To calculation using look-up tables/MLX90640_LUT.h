/**
 * @copyright (C) 2017 Melexis N.V.
 *
 */
#ifndef _MLX640_LUT_H_
#define _MLX640_LUT_H_

    int MLX90640_GenerateLUTs(float tMin, float tMax, float tStep, paramsMLX90640 *params, float *lut1, float *lut2);
    void MLX90640_LookUpTo(uint16_t *frameData, const paramsMLX90640 *params, float emissivity, float tr, float tMin, float tStep, uint16_t lutLines, float *lut1, float *lut2, float *result);
    
#endif
