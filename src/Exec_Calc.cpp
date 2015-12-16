#include "Exec_Calc.h"
#include "CpptrajStdio.h"
#include "RPNcalc.h"

void Exec_Calc::Help() const {
  mprintf("\t<expression>\n"
          "  Evaluate the given mathematical expression.\n");
}

Exec::RetType Exec_Calc::Execute(CpptrajState& State, ArgList& argIn) {
  RPNcalc calc;
  calc.SetDebug( State.Debug() );
  // Do NOT include command in expression.
  if (calc.ProcessExpression( argIn.ArgString().substr(argIn[0].size()) ))
    return CpptrajState::ERR;
  if (calc.Evaluate( State.DSL() )) return CpptrajState::ERR;
  return CpptrajState::OK;
}
