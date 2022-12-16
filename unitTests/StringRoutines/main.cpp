// Unit test for StringRoutines 
#include <cstdio>
#include "StringRoutines.h"

static const int Err(const char* msg) {
  fprintf(stderr, "Error: %s\n", msg);
  return 1;
}

int main() {
  // Test a valid mask
  if (!StrIsMask(":1-12&!@H="))
    return Err("StrIsMask failed for valid mask.");
  // Test an invalid mask
  if (StrIsMask("test@"))
    return Err("StrIsMask failed for invalid mask.");

  return 0;
}
