#include "Exec_CreateSet.h"
#include "CpptrajStdio.h"
#include "RPNcalc.h"

// Exec_CreateSet::Help()
void Exec_CreateSet::Help() const
{
  mprintf("\t<expression> [xmin <min>] xstep <step> nx <nxvals>\n");
}

// Exec_CreateSet::Execute()
Exec::RetType Exec_CreateSet::Execute(CpptrajState& State, ArgList& argIn)
{
  // Get keywords
  double xmin  = argIn.getKeyDouble("xmin", 0.0);
  double xstep = argIn.getKeyDouble("xstep", 0.0);
  if (xstep <= 0.0) {
    mprinterr("Error: 'xstep' must be specified and not be 0.0\n");
    return CpptrajState::ERR;
  }
  int nx = argIn.getKeyInt("nx", 0);
  if (nx < 1) {
    mprinterr("Error: 'nx' must be specified and > 0\n");
    return CpptrajState::ERR;
  }
  // Get equation
  std::string equation = argIn.GetStringNext();
  if (equation.empty()) {
    mprinterr("Error: Must specify an equation.\n");
    return CpptrajState::ERR;
  }
  RPNcalc Calc;
  Calc.SetDebug( State.Debug() );
  if (Calc.ProcessExpression( equation )) return CpptrajState::ERR;
  // Equation must have an assignment
  if ( Calc.AssignStatus() != RPNcalc::YES_ASSIGN ) {
    mprinterr("Error: No assignment '=' in equation '%s'\n", equation.c_str());
    return CpptrajState::ERR;
  }
  std::string dsoutName = Calc.FirstTokenName();
  if (dsoutName.empty()) {
    mprinterr("Error: Invalid assignment in equation '%s'\n", equation.c_str());
    return CpptrajState::ERR;
  }
  // Create X array
  DataSet* xset = State.DSL().AddSet( DataSet::DOUBLE, MetaData("X") );
  if (xset == 0) return CpptrajState::ERR;
  xset->Allocate( DataSet::SizeArray(1, nx) );
  double xval = xmin;
  for (unsigned int i = 0; i < (unsigned int)nx; i++) {
    xset->Add(i, &xval);
    xval += xstep;
  }
  // Evaluate
  if (Calc.Evaluate( State.DSL() )) return CpptrajState::ERR;
  // Get the resulting data set
  DataSet* dsout = State.DSL().GetDataSet( dsoutName );
  if (dsout == 0) {
    mprinterr("Error: Set '%s' not found.\n", dsoutName.c_str());
    return CpptrajState::ERR;
  }
  // Remove X array
  State.DSL().RemoveSet( xset );
  // Set output dimension
  dsout->SetDim( Dimension::X, Dimension(xmin, xstep, "X") );
  
  return CpptrajState::OK;
}
