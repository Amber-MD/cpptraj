#include "AssociatedData.h"
#include "CpptrajStdio.h"
#include "ArgList.h"

const char* AssociatedData_NOE::HelpText = 
  "[bound <lower> bound <upper>] [rexp <expected>] [noe_strong] [noe_medium] [noe_weak]";

int AssociatedData_NOE::NOE_Args(ArgList& argIn) {
  l_bound_ = argIn.getKeyDouble("bound", 0.0);
  u_bound_ = argIn.getKeyDouble("bound", 0.0);
  rexp_ = argIn.getKeyDouble("rexp", -1.0);
  if (argIn.hasKey("noe_weak")) {
    l_bound_ = 3.5;
    u_bound_ = 5.0;
  } else if (argIn.hasKey("noe_medium")) {
    l_bound_ = 2.9;
    u_bound_ = 3.5;
  } else if (argIn.hasKey("noe_strong")) {
    l_bound_ = 1.8;
    u_bound_ = 2.9;
  }
  if (u_bound_ <= l_bound_) {
    mprinterr("Error: noe lower bound (%g) must be less than upper bound (%g).\n",
              l_bound_, u_bound_);
    return 1;
  }
  return 0;
}
