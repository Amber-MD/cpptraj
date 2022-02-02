// Unit test for NameType class
#include <cstdio>
#include <iostream>
#include <cmath>

#include <string>
#include <algorithm>
#include "GistEntropyUtils.h"

static const int Err(const char* msg) {
  fprintf(stderr, "Error: %s\n", msg);
  return 1;
}
  
template<typename T>
bool close(T a, T b)
{
  T atol = 1e-8;
  T rtol = 1e-8;
  return fabs(a - b) < (atol + fabs(b) * rtol);
}

int main() {
  const double vals[] = { 0, 1, 3 };
  double spacing = 2;
  if (interpolate(1, vals, 3, spacing) != 0.5) { return Err("Wrong interpolation at x=1."); }
  if (interpolate(0, vals, 3, spacing) != 0.) { return Err("Wrong interpolation at x=0."); }
  if (interpolate(3, vals, 3, spacing) != 2.) { return Err("Wrong interpolation at x=3."); }
  // different spacing
  if (!close(interpolate(0.3, vals, 3, 0.2), 2.)) { return Err("Wrong interpolation at x=0.3 with non-integer spacing."); }
  // extrapolate linearly to the left
  if (interpolate(-4, vals, 3, spacing) != -2) { return Err("Wrong interpolation at x=-4."); }
  // extrapolate linearly to the right
  if (interpolate(8, vals, 3, spacing) != 7) { return Err("Wrong interpolation at x=8."); }
  

  if (sixVolumeCorrFactor(0) != 1.) { return Err("Wrong six corr factor at x=0."); }
  if (!close(sixVolumeCorrFactor(4), 1.785045060669792516)) { return Err("Wrong six corr factor at x=4."); }
  return 0;
}
