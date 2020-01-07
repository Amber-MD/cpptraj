#ifndef INC_FORLOOP_INTEGER_H
#define INC_FORLOOP_INTEGER_H
#include "ForLoop.h"
/// Basic for loops over integer value
class ForLoop_integer : public ForLoop {
  public:
    ForLoop_integer();

    static void helpText();

    int SetupFor(CpptrajState&, ArgList&);
    int BeginFor(DataSetList const&);
    bool EndFor(DataSetList&);

    enum OpType { INCREMENT=0, DECREMENT, LESS_THAN, GREATER_THAN, LT_EQUALS, GT_EQUALS, NO_OP };
  private:
    static const char* OpStr_[NO_OP];

    int calcNumIterations() const;

    std::string startVarName_;   ///< Variable containing initial value
    std::string endVarName_;     ///< Variable containing the end value
    OpType endOp_;               ///< (INTEGER only) end operator
    OpType incOp_;               ///< (INTEGER only) increment operator
    int start_;                  ///< (INTEGER only) initial value
    int end_;                    ///< (INTEGER only) end value
    int inc_;                    ///< (INTEGER only) increment value
    int currentVal_;             ///< (INTEGER only) current value
    bool endArgPresent_;         ///< True if an end argument was specified
};
#endif
