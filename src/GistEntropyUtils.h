#ifndef GIST_ENTROPY_UTILS_H
#define GIST_ENTROPY_UTILS_H

#include <cstdlib>

const double SIX_CORR_SPACING = 0.01;

/**
 * Linear interpolation with extrapolation for out-of-bounds values.
 * 
 * Given discrete function f(x) with *values* at positions starting from zero with given *spacing*,
 * return the linearly interpolated value at *x*.
 * If x is out of bounds, extrapolates linearly using the two first (or last) values.
 * 
 * @param  x       position to interpolate values at
 * @param  values  discrete function values
 * @param  spacing spacing of x positions
 */
double interpolate(double x, const double* values, size_t array_len, double spacing);

// Interpolate SIX_CORR at a given position
double sixVolumeCorrFactor(double NNs);

#endif