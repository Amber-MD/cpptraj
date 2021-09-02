#include <cstdio>
#include <cmath>
#include "ArgList.h"

static const int Err(const char* msg) {
  fprintf(stderr, "Error: %s\n", msg);
  return 1;
}

const double SMALL        = 0.00000000000001;

bool f_not_equals(double d1, double d2) {
  return (fabs(d1 - d2) > SMALL);
}

int main() {
  //          7       6 6    4               1    1 2  3    3       5
  ArgList a1("This is a test 3.14159 arglist with 4 45 darg 123.323 :12@H=");

  int iarg = a1.getKeyInt("with", -1);
  if (iarg != 4) return Err("getKeyInt() failed.");

  iarg = a1.getNextInteger(-1);
  if (iarg != 45) return Err("getNextInteger() failed.");

  double darg = a1.getKeyDouble("darg", 0);
  if (f_not_equals(darg, 123.323)) return Err("getKeyDouble() failed.");

  darg = a1.getNextDouble(0);
  if (f_not_equals(darg, 3.14159)) return Err("getNextDouble() failed.");

  std::string sarg = a1.GetMaskNext();
  if (sarg != ":12@H=") return Err("GetMaskNext() failed.");

  sarg = a1.GetStringKey("a");
  if (sarg != "test") return Err("GetStringKey() failed.");

  sarg = a1.GetStringNext();
  if (sarg != "This") return Err("GetStringKey() failed.");

  ArgList r1 = a1.RemainingArgs();
  if (a1.NremainingArgs() > 0) return Err("RemainingArgs() failed (# remaining args).");
  if (r1.Nargs() != 2) return Err("RemainingArgs() failed (# new args).");

  r1.MarkArg(0);
  sarg = r1.GetStringNext();
  if (sarg != "arglist") return Err("MarkArg() failed.");

  return 0;
}
