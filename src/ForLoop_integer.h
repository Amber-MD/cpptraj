#ifndef INC_FORLOOP_INTEGER_H
#define INC_FORLOOP_INTEGER_H
#include "ForLoop.h"
/// Basic for loops over integer value
class ForLoop_integer : public ForLoop {
  public:
    ForLoop_integer();

    int SetupFor(CpptrajState&, std::string const&, ArgList&);
    int BeginFor();
    bool EndFor(VariableArray&);

    enum OpType { INCREMENT=0, DECREMENT, LESS_THAN, GREATER_THAN, NO_OP };
  private:
    static const char* OpStr_[NO_OP];

    OpType endOp_;               ///< (INTEGER only) end operator
    OpType incOp_;               ///< (INTEGER only) increment operator
    int start_;                  ///< (INTEGER only) initial value
    int end_;                    ///< (INTEGER only) end value
    int inc_;                    ///< (INTEGER only) increment value
    int currentVal_;             ///< (INTEGER only) current value
};
#endif
