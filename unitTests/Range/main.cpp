// Unit test for Range class
#include <cstdio>
#include "Range.h"

static const int Err(const char* msg) {
  fprintf(stderr, "Error: %s\n", msg);
  return 1;
}

int main() {
  Range myRange;
  if (myRange.SetRange("3-6,9,12-14")) return Err("SetRange failed.");
  myRange.PrintRange("myRange + 1", 1);
  printf("\n");
  printf("Range arg: '%s'\n", myRange.RangeArg());
  Range::const_iterator it = myRange.begin();
  for (int i = 3; i < 7; i++, ++it) {
    if (*it != i) return Err("Iterator range 3-6 does not match.\n");
  }
  if (*it != 9) return Err("Iterator 9 does not match.\n");
  ++it;
  for (int i = 12; i < 15; i++, ++it) {
    if (*it != i) return Err("Iterator range 12-14 does not match.\n");
  }

  return 0;
}
