#ifndef INC_DISPATCHOBJECT_H
#define INC_DISPATCHOBJECT_H
/// Abstract base class that all dispatchable objects will inherit.
class DispatchObject {
  public:
    enum DispatchType { NONE=0, PARM, COORD, ACTION, ANALYSIS, GENERAL, DEPRECATED };
    // Function pointers
    typedef DispatchObject* (*DispatchAllocatorType)();
    typedef void (*DispatchHelpType)();
    //typedef bool (*DispatchCmdType)(ArgList const&);
    typedef const char* DispatchCmdType;
    // Struct that describes how a dispatchable object is called.
    struct Token {
      DispatchType Type;
      DispatchCmdType Cmd;
      DispatchAllocatorType Alloc;
      DispatchHelpType Help;
      int Idx; // TODO: Remove
    };
    typedef const Token* TokenPtr;
};
#endif
