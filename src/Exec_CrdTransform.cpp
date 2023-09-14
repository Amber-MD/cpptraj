#include "Exec_CrdTransform.h"
#include "CpptrajStdio.h"
#include "Action_Align.h"

// Exec_CrdTransform::Help()
void Exec_CrdTransform::Help() const
{

}

/** Transform coordinates by RMS-fitting to an average structure, calculating
  * a new average, then RMS-fitting to that average and so on until a
  * tolerance is reached. Essentially the procedure described by 
  * Klem et al. J. Chem. Theory Comput. 2022, 18, 3218−3230.
  */
int Exec_CrdTransform::iterativeRmsRefinement(CpptrajState& State,
                                              AtomMask const& maskIn,
                                              bool useMass,
                                              double tolIn,
                                              DataSet_Coords* crdIn,
                                              DataSet_Coords* crdOut)
const
{
  mprintf("\tRMS iterative refinement.\n");
  mprintf("\tInput coords: %s\n", crdIn->legend());
  mprintf("\tOutput coords: %s\n", crdIn->legend());
  mprintf("\tAtom mask: %s\n", maskIn.MaskString());
  mprintf("\tRMS Tolerance: %g Ang.\n", tolIn);
  if (useMass)
    mprintf("\tMass-weighting on.\n");
  else
    mprintf("\tMass-weighting off.\n");
  // Do the initial fit to the first frame.
  Frame frmIn = crdIn->AllocateFrame();
  crdIn->GetFrame(0, frmIn);
  Frame selectedRef;
  selectedRef.SetupFrameFromMask( maskIn, crdIn->Top().Atoms() );
  selectedRef.SetCoordinates( frmIn, maskIn );
  // Ensure reference is centered on the origin
  Vec3 refTrans = selectedRef.CenterOnOrigin( useMass );
  // Set up frame for selected incoming atoms
  Frame selectedTgt = selectedRef;
  // Set up frame to hold average 
  Frame avgFrm = selectedTgt;

  double currentTol = tolIn + 9999.0;

  while (currentTol > tolIn) {
    avgFrm.ZeroCoords();
    Vec3 tgtTrans(0.0);
    Matrix_3x3 rot(0.0);
    for (unsigned int idx = 0; idx != crdIn->Size(); idx++) {
      crdIn->GetFrame(idx, frmIn);
      selectedTgt.SetCoordinates( frmIn, maskIn );
      selectedTgt.RMSD_CenteredRef( selectedRef, rot, tgtTrans, useMass );
      frmIn.Trans_Rot_Trans(tgtTrans, rot, refTrans);
      crdOut->SetCRD(idx, frmIn);
      avgFrm.AddByMask( frmIn, maskIn );
    }
    avgFrm.Divide( (double)crdIn->Size() );
    // Calc RMS of current average to current reference
    currentTol = avgFrm.RMSD_CenteredRef( selectedRef, rot, tgtTrans, useMass );
    // Fit the current average TODO is this necessary?
    avgFrm.Trans_Rot_Trans(tgtTrans, rot, refTrans);
    // Set current average to be new reference
    selectedRef = avgFrm;
  }

  return 0;
}
  

// Exec_CrdTransform::Execute()
Exec::RetType Exec_CrdTransform::Execute(CpptrajState& State, ArgList& argIn)
{
  AtomMask mask( argIn.GetStringKey("mask") );
  bool useMass = argIn.hasKey("mass");
  double rmsTol = argIn.getKeyDouble("rmstol", 0.0001);
  
  // Get COORDS set
  std::string setname = argIn.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: %s: Specify COORDS dataset name.\n", argIn.Command());
    return CpptrajState::ERR;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)State.DSL().FindSetOfGroup( setname, DataSet::COORDINATES );
  if (CRD == 0) {
    mprinterr("Error: %s: No COORDS set with name %s found.\n", argIn.Command(), setname.c_str());
    return CpptrajState::ERR;
  }
  mprintf("\tUsing set '%s'\n", CRD->legend());
  if (CRD->Size() < 1) {
    mprinterr("Error: Set '%s' has no frames.\n", CRD->legend());
    return CpptrajState::ERR;
  }
  if (CRD->Type() == DataSet::TRAJ) {
    mprinterr("Error: TRAJ sets not yet supported.\n"); // FIXME
    return CpptrajState::ERR;
  }

  // Set up mask
  if (CRD->Top().SetupIntegerMask( mask )) {
    mprinterr("Error: Could not set up mask.\n");
    return CpptrajState::ERR;
  }
  mask.MaskInfo();

  // RMS iterative refinement
  int err = iterativeRmsRefinement(State, mask, useMass, rmsTol, CRD, CRD);
  if (err != 0) return CpptrajState::ERR;

  return CpptrajState::OK;
}
