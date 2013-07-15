#ifndef INC_COMMANDLIST_H
#define INC_COMMANDLIST_H
#include "CpptrajState.h"
#include "ArgList.h"
#include "DispatchObject.h"
//enum CmdReturnType { C_OK = 0, C_ERR, C_QUIT, C_INTERACTIVE };
class CommandList {
  public:
    /// Unique command IDs TODO: Are these necessary?
    enum CommandID { 
      NO_ID = 0, LIST = 0, HELP, QUIT, RUN, DEBUG, NOPROG, NOEXITERR, SYSTEM,
      ACTIVEREF, READDATA, CREATE, PRECISION, DATAFILE, SELECT, SELECTDS,
      READINPUT, RUN_ANALYSIS, WRITEDATA, CLEAR, LOADCRD, CRDACTION, CRDOUT,
      WRITE,
      // TRAJ
      REFERENCE, TRAJIN, TRAJOUT,
      // PARM
      LOADPARM, PARMINFO, PARMWRITE, PARMSTRIP, PARMBOX, SOLVENT, BONDINFO,
      RESINFO, MOLINFO, CHARGEINFO, SCALEDIHEDRALK
    };
    // TODO: Make below private, make commands part of CommandList class?
    /// Command categories.
    enum CommandType { NONE=0, PARM, TRAJ, ACTION, ANALYSIS, GENERAL, DEPRECATED };
    /// Function pointer to command function.
    typedef int (*CommandFxnType)(CpptrajState&, ArgList&,
                                  DispatchObject::DispatchAllocatorType, int);
    /// Function pointer to help function.
    typedef void (*CommandHelpType)();
    /// Keyword type.
    typedef const char* CommandKeywordType;
    /// Struct that describes how a command is called.
    struct Token {
      CommandType Type;                            ///< Command type
      CommandKeywordType Cmd;                      ///< Command keyword
      DispatchObject::DispatchAllocatorType Alloc; ///< Allocator (Action/Analysis only)
      CommandHelpType Help;                        ///< Help text function.
      CommandID Idx;                               ///< Command ID (all except Action/Analysis).
      CommandFxnType Fxn;                          ///< Command function.
    };
    /// Pointer to command token.
    typedef const Token* TokenPtr;

    static void List(CommandType);
    static TokenPtr SearchTokenType(CommandType, ArgList const& argIn);
    static TokenPtr SearchToken(ArgList&);
  private:
    
    static const char* CommandTitle[];
    /// Master list of commands.
    static const Token Commands[];

    int debug_;
};
#endif
