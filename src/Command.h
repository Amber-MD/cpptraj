#ifndef INC_COMMAND_H
#define INC_COMMAND_H
#include "DispatchObject.h"
#include "ArgList.h"
/// This is a static class that determines how commands are handled.
/** To add a new Action/Analysis command, add the appropriate '#include'
  * to the top of Command.cpp and add an entry to Commands[] (search for 
  * INC_ACTION/INC_ANALYSIS respectively).
  */
class Command {
  public:
    enum Type { LIST = 0, HELP, QUIT, RUN, DEBUG, NOPROG, NOEXITERR,
                SYSTEM, ACTIVEREF, READDATA, CREATE, PRECISION, DATAFILE,
                SELECT, SELECTDS, READINPUT, RUN_ANALYSIS, WRITEDATA,
                CLEAR, LOADCRD, CRDACTION, CRDOUT, WRITE,
                // TRAJ
                REFERENCE, TRAJIN, TRAJOUT,
                // PARM
                LOADPARM, PARMINFO, PARMWRITE, PARMSTRIP, PARMBOX,
                SOLVENT, BONDINFO, RESINFO, MOLINFO, CHARGEINFO };
    static void List(DispatchObject::DispatchType);
    static DispatchObject::TokenPtr SearchToken(ArgList&);
    static DispatchObject::TokenPtr SearchTokenType(DispatchObject::DispatchType, ArgList&);
  private:
    static const DispatchObject::Token Commands[];
};
#endif 
