#ifndef INC_DISPATCHOBJECT_H
#define INC_DISPATCHOBJECT_H
/// Abstract base class that all dispatchable objects will inherit.
class DispatchObject {
  public:
    enum DispatchType { NONE=0, PARM, TRAJ, ACTION, ANALYSIS, GENERAL, DEPRECATED };
    // Function pointers
    typedef DispatchObject* (*DispatchAllocatorType)();
    typedef void (*DispatchHelpType)();
    //typedef bool (*DispatchCmdType)(ArgList const&);
    typedef const char* DispatchCmdType;
    /// Struct that describes how a command is called.
    struct Token {
      DispatchType Type;           ///< Command type
      DispatchCmdType Cmd;         ///< Command keyword
      DispatchAllocatorType Alloc; ///< Allocator (Action/Analysis only)
      DispatchHelpType Help;       ///< Help text function.
      int Idx;                     ///< Command index (all except Action/Analysis).
    };
    typedef const Token* TokenPtr;
};
#endif
