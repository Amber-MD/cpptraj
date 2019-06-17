#include "Action_DihedralRMS.h"
#include "CpptrajStdio.h"
#include "TorsionRoutines.h"

// Action_DihedralRMS::Help()
void Action_DihedralRMS::Help() const {
  mprintf("\t[<name>] <dihedral types> [out <file>]\n"
          "%s"
          "\t[%s]\n", ReferenceAction::Help(), DihedralSearch::newTypeArgsHelp_);
}

// Action_DihedralRMS::Init()
Action::RetType Action_DihedralRMS::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get keywords
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs);
  std::string tgtArg = actionArgs.GetStringKey("tgtrange");
  std::string refArg = actionArgs.GetStringKey("refrange");
  if (!tgtArg.empty()) {
    if (tgtRange_.SetRange( tgtArg )) return Action::ERR;
  }
  if (!refArg.empty()) {
    if (refRange_.SetRange( refArg )) return Action::ERR;
  } else if (!tgtArg.empty()) {
    if (refRange_.SetRange( tgtArg )) return Action::ERR;
  }
  // Reference keywords. false= no fit, false= no use mass
  if (REF_.InitRef(actionArgs, init.DSL(), false, false)) return Action::ERR;
  // Search for known dihedral keywords
  dihSearch_.SearchForArgs(actionArgs);
  // Get custom dihedral arguments: dihtype <name>:<a0>:<a1>:<a2>:<a3>[:<offset>]
  if (dihSearch_.SearchForNewTypeArgs(actionArgs)) return Action::ERR;
  // If no dihedral types yet selected, this will select all.
  dihSearch_.SearchForAll();
  // Setup DataSet(s) name
  std::string dsetname = actionArgs.GetStringNext();
  // Setup output data set
  dataOut_ = init.DSL().AddSet(DataSet::DOUBLE, dsetname, "DIHRMS");
  if (dataOut_ == 0) return Action::ERR;
  if (outfile != 0) outfile->AddDataSet( dataOut_ );
# ifdef MPI
  if (REF_.SetTrajComm( init.TrajComm() )) return Action::ERR;
# endif

  mprintf("    DIHEDRAL RMSD: Calculating dihedral RMS for dihedrals:");
  dihSearch_.PrintTypes();
  mprintf("\n");
  if (!tgtRange_.Empty())
    mprintf("\tTarget residue range: %s\n", tgtRange_.RangeArg());
  if (!refRange_.Empty())
    mprintf("\tReference residue range: %s\n", refRange_.RangeArg());
  mprintf("\tReference is %s\n", REF_.RefModeString().c_str());
  return Action::OK;
}

/** Fill the refVals_ array with reference dihedral values. */
int Action_DihedralRMS::GetRefDihedrals(Topology const& top, Frame const& frm) {
  Range actualRange = GetActualRange(top, refRange_);
  if (actualRange.Empty()) {
    mprinterr("Error: No residues in reference topology %s\n", top.c_str());
    return 1;
  }
  // Search for specified dihedrals in each residue in the range
  if (dihSearch_.FindDihedrals(top, actualRange)) {
    mprinterr("Error: No dihedrals found in reference topology %s\n", top.c_str());
    return 1;
  }
  mprintf("DEBUG: %i reference dihedrals\n", dihSearch_.Ndihedrals());
  // Calculate dihedrals
  refVals_.clear();
  refVals_.reserve( dihSearch_.Ndihedrals() );
  for (DihedralSearch::mask_it dih = dihSearch_.begin();
                               dih != dihSearch_.end(); ++dih)
  {
    double torsion = Torsion( frm.XYZ(dih->A0()),
                              frm.XYZ(dih->A1()),
                              frm.XYZ(dih->A2()),
                              frm.XYZ(dih->A3()) );
    refVals_.push_back( torsion );
  }
  return 0;
}

/** \return Range with actual residues to search for dihedrals in. If the
  *         given range is not set, will use all solute residues.
  */
Range Action_DihedralRMS::GetActualRange(Topology const& top, Range const& resRange) const {
  Range actualRange;
  // If range is empty (i.e. no resrange arg given) look through all 
  // solute residues.
  if (resRange.Empty())
    actualRange = top.SoluteResidues();
  else {
    // If user range specified, create new range shifted by -1 since internal
    // resnums start from 0.
    actualRange = resRange;
    actualRange.ShiftBy(-1);
  }
  return actualRange;
}

// Action_DihedralRMS::Setup()
Action::RetType Action_DihedralRMS::Setup(ActionSetup& setup)
{
  // Perform any reference setup. -1 means do not check number of selected
  // atoms, which is not necessary because we are looking for dihedrals,
  // not atoms.
  if (REF_.SetupRef(setup.Top(), -1))
    return Action::SKIP;
  // Fill the reference values if needed
  if (refVals_.empty() && REF_.RefMode() == ReferenceAction::FRAME) {
    // Sanity check
    Topology* refTop = REF_.RefCrdTopPtr();
    if (refTop == 0) {
      mprinterr("Internal Error: Reference topology is null.\n");
      return Action::ERR;
    }
    if (GetRefDihedrals( *refTop, REF_.CurrentReference() )) return Action::ERR;
    // DEBUG
    mprintf("DEBUG: Reference dihedral values (radians):\n");
    for (Darray::const_iterator it = refVals_.begin(); it != refVals_.end(); ++it)
      mprintf("\t%8li %12.4f\n", it-refVals_.begin(), *it);
  }
  // Get the target range
  Range actualRange = GetActualRange(setup.Top(), tgtRange_);
  if (actualRange.Empty()) {
    mprinterr("Error: No residues in topology %s\n", setup.Top().c_str());
    return Action::ERR;
  }
  // Search for specified dihedrals in each residue in the range
  if (dihSearch_.FindDihedrals(setup.Top(), actualRange)) {
    mprinterr("Error: No dihedrals found in topology %s\n", setup.Top().c_str());
    return Action::ERR;
  }
  mprintf("\tSelected %i dihedrals\n", dihSearch_.Ndihedrals());

  return Action::OK;
}

// Action_DihedralRMS::DoAction()
Action::RetType Action_DihedralRMS::DoAction(int frameNum, ActionFrame& frm)
{

  return Action::ERR;
}
