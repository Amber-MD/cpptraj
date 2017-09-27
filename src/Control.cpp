#include "Control.h"
#include "CpptrajStdio.h"

void Control_For::Help() const {
  mprintf("\tTODO\n");
}

int Control_For::SetupControl(ArgList& argIn) {
  // for <var> in <mask> do
  std::string in_arg = argIn.GetStringKey("in");
  if (in_arg.empty()) {
    mprinterr("Error: 'for': expected 'in <expression>'\n");
    return 1;
  }
  if (!argIn.hasKey("do")) {
    mprinterr("Error: 'for': missing 'do'\n");
    return 1;
  }
  varname_ = argIn.GetStringNext();
  if (varname_.empty()) {
    mprinterr("Error: 'for': missing variable name.\n");
    return 1;
  }
  mprintf("DEBUG: for %s in %s do\n", varname_.c_str(), in_arg.c_str());
  return 0;
}
