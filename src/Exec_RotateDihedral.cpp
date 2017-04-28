#include "Exec_RotateDihedral.h"
#include "CpptrajStdio.h"
#include "DihedralSearch.h"
#include "Constants.h"
#include "TorsionRoutines.h"

void Exec_RotateDihedral::Help() const {
  mprintf("\tcrdset <COORDS set> [frame <#>] [name <output set name>]\n"
          "\t{value <value> | increment <increment>}\n"
          "\t{ <mask1> <mask2> <mask3> <mask4> |\n"
          "\t  res <#> type <dih type> }\n"
          "\t<dih type> =");
  DihedralSearch::ListKnownTypes();
  mprintf("  Rotate specified dihedral to specified value or by given increment.\n");
}

static const char* ModeStr[] = { "value", "increment" };

static inline Exec::RetType MaskError(AtomMask const& mask) {
  mprinterr("Error: Mask '%s' selects %i atoms, expected 1.\n",
            mask.MaskString(), mask.Nselected());
  return CpptrajState::ERR;
}

// Exec_RotateDihedral::Execute()
Exec::RetType Exec_RotateDihedral::Execute(CpptrajState& State, ArgList& argIn) {
  // Get input COORDS set
  std::string setname = argIn.GetStringKey("crdset");
  if (setname.empty()) {
    mprinterr("Error: Specify COORDS dataset name with 'crdset'.\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)State.DSL().FindCoordsSet( setname );
  if (CRD == 0) {
    mprinterr("Error: Could not find COORDS set '%s'\n", setname.c_str());
    return CpptrajState::ERR;
  }
  if (CRD->Size() < 1) {
    mprinterr("Error: COORDS set is empty.\n");
    return CpptrajState::ERR;
  }
  int frame = argIn.getKeyInt("frame", 0);
  if (frame < 0 || frame >= (int)CRD->Size()) {
    mprinterr("Error: Specified frame %i is out of range.\n", frame+1);
    return CpptrajState::ERR;
  }
  mprintf("    ROTATEDIHEDRAL: Using COORDS '%s', frame %i\n", CRD->legend(), frame+1);
  // Get target frame
  Frame FRM = CRD->AllocateFrame();
  CRD->GetFrame(frame, FRM);
  // Save as reference
  Frame refFrame = FRM;

  // Create output COORDS set if necessary
  DataSet_Coords* OUT = 0;
  int outframe = 0;
  std::string outname = argIn.GetStringKey("name");
  if (outname.empty()) {
    // This will not work for TRAJ data sets
    if (CRD->Type() == DataSet::TRAJ) {
      mprinterr("Error: Using TRAJ as input set requires use of 'name' keyword for output.\n");
      return CpptrajState::ERR;
    }
    OUT = CRD;
    outframe = frame;
  } else {
    // Create new output set with 1 empty frame.
    OUT = (DataSet_Coords*)State.DSL().AddSet( DataSet::COORDS, outname );
    if (OUT == 0) return CpptrajState::ERR;
    OUT->Allocate( DataSet::SizeArray(1, 1) );
    OUT->CoordsSetup( CRD->Top(), CRD->CoordsInfo() );
    OUT->AddFrame( CRD->AllocateFrame() );
    mprintf("\tOutput to set '%s'\n", OUT->legend());
  }

  // Determine whether we are setting or incrementing.
  enum ModeType { SET = 0, INCREMENT };
  ModeType mode = SET;
  if (argIn.Contains("value"))
    mode = SET;
  else if (argIn.Contains("increment"))
    mode = INCREMENT;
  else {
    mprinterr("Error: Specify 'value <value>' or 'increment <increment>'\n");
    return CpptrajState::ERR;
  }
  double value = argIn.getKeyDouble(ModeStr[mode], 0.0);
  switch (mode) {
    case SET: mprintf("\tDihedral will be set to %g degrees.\n", value); break;
    case INCREMENT: mprintf("\tDihedral will be incremented by %g degrees.\n", value); break;
  }
  // Convert to radians
  value *= Constants::DEGRAD;

  // Select dihedral atoms
  int A1, A2, A3, A4;
  if (argIn.Contains("type")) {
    // By type
    ArgList typeArg = argIn.GetStringKey("type");
    if (typeArg.empty()) {
      mprinterr("Error: No dihedral type specified after 'type'\n");
      return CpptrajState::ERR;
    }
    DihedralSearch dihSearch;
    dihSearch.SearchForArgs( typeArg );
    if (dihSearch.NoDihedralTokens()) {
      mprinterr("Error: Specified dihedral type not recognized.\n");
      return CpptrajState::ERR;
    }
    // Get residue
    int res = argIn.getKeyInt("res", -1);
    if (res <= 0) {
      mprinterr("Error: If 'type' specified 'res' must be specified and > 0.\n");
      return CpptrajState::ERR;
    }
    // Search for dihedrals. User residue #s start from 1.
    if (dihSearch.FindDihedrals(CRD->Top(), Range(res-1)))
      return CpptrajState::ERR;
    DihedralSearch::mask_it dih = dihSearch.begin();
    A1 = dih->A0();
    A2 = dih->A1();
    A3 = dih->A2();
    A4 = dih->A3();
  } else {
    // By masks
    AtomMask m1( argIn.GetMaskNext() );
    AtomMask m2( argIn.GetMaskNext() );
    AtomMask m3( argIn.GetMaskNext() );
    AtomMask m4( argIn.GetMaskNext() );
    if (CRD->Top().SetupIntegerMask( m1 )) return CpptrajState::ERR;
    if (CRD->Top().SetupIntegerMask( m2 )) return CpptrajState::ERR;
    if (CRD->Top().SetupIntegerMask( m3 )) return CpptrajState::ERR;
    if (CRD->Top().SetupIntegerMask( m4 )) return CpptrajState::ERR;
    if (m1.Nselected() != 1) return MaskError( m1 );
    if (m2.Nselected() != 1) return MaskError( m2 );
    if (m3.Nselected() != 1) return MaskError( m3 );
    if (m4.Nselected() != 1) return MaskError( m4 );
    A1 = m1[0];
    A2 = m2[0];
    A3 = m3[0];
    A4 = m4[0];
  }
  mprintf("\tRotating dihedral defined by atoms '%s'-'%s'-'%s'-'%s'\n",
          CRD->Top().AtomMaskName(A1).c_str(),
          CRD->Top().AtomMaskName(A2).c_str(),
          CRD->Top().AtomMaskName(A3).c_str(),
          CRD->Top().AtomMaskName(A4).c_str());
  // Set mask of atoms that will move during dihedral rotation
  AtomMask Rmask = DihedralSearch::MovingAtoms(CRD->Top(), A2, A3);
  // Calculate current value of dihedral
  double torsion = Torsion( FRM.XYZ(A1), FRM.XYZ(A2), FRM.XYZ(A3), FRM.XYZ(A4) );
  // Calculate delta needed to get to target value.
  double delta;
  switch (mode) {
    case SET:       delta = value - torsion; break;
    case INCREMENT: delta = value; break;
  }
  mprintf("\tOriginal torsion is %g, rotating by %g degrees.\n",
          torsion*Constants::RADDEG, delta*Constants::RADDEG);
  // Set axis of rotation
  Vec3 axisOfRotation = FRM.SetAxisOfRotation( A2, A3 );
  // Calculate rotation matrix for delta.
  Matrix_3x3 rotationMatrix;
  rotationMatrix.CalcRotationMatrix(axisOfRotation, delta);
  // Rotate around axis
  FRM.Rotate(rotationMatrix, Rmask);
  // RMS-fit the non-moving part of the coords back on original
  AtomMask refMask = Rmask;
  refMask.InvertMask();
  FRM.Align( refFrame, refMask );
  // Update coords
  OUT->SetCRD( outframe, FRM );

  return CpptrajState::OK;
}
