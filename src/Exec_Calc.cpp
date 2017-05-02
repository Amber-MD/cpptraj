#include "Exec_Calc.h"
#include "CpptrajStdio.h"
#include "RPNcalc.h"

void Exec_Calc::Help() const {
  mprintf("\t<expression>\n\t[prec <width>.<precision>] [format {double|general|scientific}]\n"
          "  Evaluate the given mathematical expression.\n");
}

Exec::RetType Exec_Calc::Execute(CpptrajState& State, ArgList& argIn) {
  RPNcalc calc;
  calc.SetDebug( State.Debug() );
  if (calc.ProcessOptions(argIn)) return CpptrajState::ERR;
  // Do NOT include command in expression.
  if (calc.ProcessExpression( argIn.ArgString() ))
    return CpptrajState::ERR;
  if (calc.Evaluate( State.DSL() )) return CpptrajState::ERR;
  return CpptrajState::OK;
}
