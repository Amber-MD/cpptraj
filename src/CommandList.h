#ifndef INC_COMMANDLIST_H
#define INC_COMMANDLIST_H
#include "CpptrajState.h"
#include "ArgList.h"
#include "DispatchObject.h"
class CommandList {
  public:
    enum RetType { C_OK = 0, C_ERR, C_QUIT, C_INTERACTIVE };
    // TODO: Make below private, make commands part of CommandList class?
    /// Command categories.
    enum CommandType { NONE=0, PARM, TRAJ, ACTION, ANALYSIS, GENERAL, DEPRECATED };
    /// Shorthand for DispatchAllocatorType
    typedef DispatchObject::DispatchAllocatorType AllocType;
    /// Function pointer to command function.
    typedef RetType (*CommandFxnType)(CpptrajState&, ArgList&, AllocType);
    /// Function pointer to help function.
    typedef void (*CommandHelpType)();
    /// Keyword type.
    typedef const char* CommandKeywordType;
    /// Struct that describes how a command is called.
    struct Token {
      CommandType Type;       ///< Command type
      CommandKeywordType Cmd; ///< Command keyword
      AllocType Alloc;        ///< Allocator (Action/Analysis only)
      CommandHelpType Help;   ///< Help text function.
      CommandFxnType Fxn;     ///< Command function.
    };
    /// Pointer to command token.
    typedef const Token* TokenPtr;

    static void ListCommands(CommandType);
    static TokenPtr SearchTokenType(CommandType, ArgList const& argIn);
    static TokenPtr SearchToken(ArgList&);
    static RetType Dispatch(CpptrajState&, std::string const&);
    static RetType ProcessInput(CpptrajState&, std::string const&);
  private:
    
    static const char* CommandTitle[];
    /// Master list of commands.
    static const Token Commands[];

    int debug_;
};
#endif
