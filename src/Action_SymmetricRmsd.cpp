#include "Hungarian.h" 
#include "Action_SymmetricRmsd.h"
#include "CpptrajStdio.h"
#include "AtomMap.h"
#include "DistRoutines.h"

// CONSTRUCTOR
Action_SymmetricRmsd::Action_SymmetricRmsd() :
  rmsd_(0), action_return_(Action::OK), remap_(false) {}

void Action_SymmetricRmsd::Help() const {
  mprintf("\t[<name>] <mask> [<refmask>] [out <filename>] [nofit] [mass] [remap]\n"
          "\t[ first | %s |\n"
          "\t  reftraj <trajname> [parm <parmname> | parmindex <#>] ]\n"
          "  Perform symmetry-corrected RMSD calculation. If 'remap' is specified\n"
          "  frames will be modified for symmetry as well.\n", DataSetList::RefArgs);
}

// Action_SymmetricRmsd::Init()
Action::RetType Action_SymmetricRmsd::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Check for keywords
  bool fit = !actionArgs.hasKey("nofit");
  bool useMass = actionArgs.hasKey("mass");
  DataFile* outfile = init.DFL().AddDataFile(actionArgs.GetStringKey("out"), actionArgs);
  remap_ = actionArgs.hasKey("remap");
  // Reference keywords
  bool previous = actionArgs.hasKey("previous");
  bool first = actionArgs.hasKey("first");
  ReferenceFrame REF = init.DSL().GetReferenceFrame( actionArgs );
  std::string reftrajname = actionArgs.GetStringKey("reftraj");
  Topology* RefParm = init.DSL().GetTopology( actionArgs );
  // Get the RMS mask string for target
  std::string tMaskExpr = actionArgs.GetMaskNext();
  if (tgtMask_.SetMaskString( tMaskExpr )) return Action::ERR;
  // Initialize Symmetric RMSD calc.
  if (SRMSD_.InitSymmRMSD( fit, useMass, debugIn )) return Action::ERR;
  // Initialize reference
  std::string rMaskExpr = actionArgs.GetMaskNext();
  if (rMaskExpr.empty())
    rMaskExpr = tMaskExpr;
  if (REF_.InitRef(previous, first, useMass, fit, reftrajname, REF, RefParm,
                   rMaskExpr, actionArgs, "symmrmsd"))
    return Action::ERR;
  // Set up the RMSD data set.
  MetaData md(actionArgs.GetStringNext(), MetaData::M_RMS); 
  rmsd_ = init.DSL().AddSet(DataSet::DOUBLE, md, "RMSD");
  if (rmsd_==0) return Action::ERR;
  // Add dataset to data file list
  if (outfile != 0) outfile->AddDataSet( rmsd_ );
  if (remap_ || SRMSD_.Fit())
    action_return_ = Action::MODIFY_COORDS;
  else
    action_return_ = Action::OK;
  
  mprintf("    SYMMRMSD: (%s), reference is %s", tgtMask_.MaskString(),
          REF_.RefModeString());
  if (!SRMSD_.Fit())
    mprintf(", no fitting");
  else
    mprintf(", with fitting");
  if (SRMSD_.UseMass())
    mprintf(", mass-weighted");
  mprintf(".\n");
  if (remap_) mprintf("\tAtoms will be re-mapped for symmetry.\n");
  return Action::OK;
}

// Action_SymmetricRmsd::Setup()
Action::RetType Action_SymmetricRmsd::Setup(ActionSetup& setup) {
  // Setup target mask.
  if (setup.Top().SetupIntegerMask( tgtMask_ )) return Action::ERR;
  tgtMask_.MaskInfo();
  if (tgtMask_.None()) {
    mprintf("Warning: No atoms selected by mask '%s'\n", tgtMask_.MaskString());
    return Action::SKIP;
  }
  // Allocate space for selected atoms in target frame. This will also
  // put the correct masses in based on the mask.
  selectedTgt_.SetupFrameFromMask(tgtMask_, setup.Top().Atoms());
  // Setup Symmetric RMSD calc (target mask, symmetric atoms etc)
  if (SRMSD_.SetupSymmRMSD( setup.Top(), tgtMask_, remap_ )) return Action::ERR;
  if (remap_) {
    // Allocate space for remapped frame; same # atoms as original frame
    remapFrame_.SetupFrameV( setup.Top().Atoms(), setup.CoordInfo() );
    targetMap_.resize( setup.Top().Natom() );
  }
  // Reference frame setup
  if (REF_.SetupRef(setup.Top(), tgtMask_.Nselected(), "symmrmsd"))
    return Action::ERR;
  return Action::OK;
}

// Action_SymmetricRmsd::DoAction()
Action::RetType Action_SymmetricRmsd::DoAction(int frameNum, ActionFrame& frm) {
  // Perform any needed reference actions
  REF_.ActionRef( frm.Frm(), SRMSD_.Fit(), SRMSD_.UseMass() );
  // Calculate symmetric RMSD
  selectedTgt_.SetCoordinates( frm.Frm(), tgtMask_ );
  double rmsdval = SRMSD_.SymmRMSD_CenteredRef( selectedTgt_, REF_.SelectedRef() );
  rmsd_->Add(frameNum, &rmsdval);
  if (remap_) {
    // Now re-map the target frame
    for (int atom = 0; atom < (int)targetMap_.size(); atom++)
      targetMap_[atom] = atom;
    SymmetricRmsdCalc::Iarray const& AMap = SRMSD_.AMap();
    for (unsigned int ref = 0; ref < AMap.size(); ++ref)
      targetMap_[ tgtMask_[ref] ] = tgtMask_[AMap[ref]];
    remapFrame_.SetCoordinatesByMap( frm.Frm(), targetMap_ );
    frm.SetFrame( &remapFrame_ );
  }
  if ( SRMSD_.Fit() )
    frm.ModifyFrm().Trans_Rot_Trans( SRMSD_.TgtTrans(), SRMSD_.RotMatrix(), REF_.RefTrans() );

  return action_return_;
}
