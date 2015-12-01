#ifndef INC_CMD_H
#define INC_CMD_H
#include "CpptrajState.h"
#include "DispatchObject.h"
/// This namespace contains the definition of a command in Cpptraj
namespace Cmd {
  /// Possible command return types.
  enum RetType { OK = 0, ERR, QUIT };
  /// Command categories.
  enum Ctype { NONE=0,  PARM,   TRAJ,   COORDS, ACTION, ANALYSIS,
               GENERAL, SYSTEM, HIDDEN, DEPRECATED };
  /// Shorthand for DispatchAllocatorType
  typedef DispatchObject::DispatchAllocatorType AllocType;
  /// Function pointer to command function.
  typedef RetType (*FxnType)(CpptrajState&, ArgList&, AllocType);
  /// Function pointer to help function.
  typedef void (*HelpType)();
  /// Keyword type.
  typedef const char* KeywordType;
  /// Struct that describes how a command is called.
  struct Token {
    Ctype Type;      ///< Command type
    KeywordType Cmd; ///< Command keyword
    AllocType Alloc; ///< Allocator (Action/Analysis only)
    HelpType Help;   ///< Help text function.
    FxnType Fxn;     ///< Command function.
  };
  /// Pointer to command token.
  typedef const Token* TokenPtr;
}
#endif
