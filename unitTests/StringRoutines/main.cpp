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
  // Test wildcard matching
  int is_match = WildcardMatch("E*.VU", "Eh5md.VU");
  if (is_match != 1)
    return Err("WildcardMatch E*.VU Eh5md.VU failed.");
  is_match = WildcardMatch("E*.VU", "Eh5md.V");
  if (is_match != 0)
    return Err("WildcardMatch E*.VU Eh5md.V failed.");
  is_match = WildcardMatch("E*.VU", "Egz.h5md.VU");
  if (is_match != 1)
    return Err("WildcardMatch E*.VU Egz.h5md.VU failed.");
  is_match = WildcardMatch("E*.VU", "Egz.h5md.V");
  if (is_match != 0)
    return Err("WildcardMatch E*.VU Egz.h5md.V failed.");

  return 0;
}
