#include "Hungarian.h" 
#include "Action_SymmetricRmsd.h"
#include "CpptrajStdio.h"
#include "AtomMap.h"
#include "DistRoutines.h"

// CONSTRUCTOR
Action_SymmetricRmsd::Action_SymmetricRmsd() : remap_(false), rmsd_(0) {}

void Action_SymmetricRmsd::Help() {
  mprintf("\t[<name>] <mask> [<refmask>] [out <filename>] [nofit] [mass] [remap]\n"
          "\t[ first | %s |\n"
          "\t  reftraj <trajname> [parm <parmname> | parmindex <#>] ]\n"
          "  Perform symmetry-corrected RMSD calculation. If 'remap' is specified\n"
          "  frames will be modified for symmetry as well.\n", FrameList::RefArgs);
}

// Action_SymmetricRmsd::Init()
Action::RetType Action_SymmetricRmsd::Init(ArgList& actionArgs, TopologyList* PFL, 
                          FrameList* FL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Check for keywords
  bool fit = !actionArgs.hasKey("nofit");
  bool useMass = actionArgs.hasKey("mass");
  DataFile* outfile = DFL->AddDataFile(actionArgs.GetStringKey("out"), actionArgs);
  remap_ = actionArgs.hasKey("remap");
  // Reference keywords
  bool previous = actionArgs.hasKey("previous");
  bool first = actionArgs.hasKey("first");
  ReferenceFrame  REF = FL->GetFrameFromArgs( actionArgs );
  std::string reftrajname = actionArgs.GetStringKey("reftraj");
  Topology* RefParm = PFL->GetParm( actionArgs );
  // Get the RMS mask string for target
  std::string tMaskExpr = actionArgs.GetMaskNext();
  // Initialize Symmetric RMSD calc.
  if (SRMSD_.InitSymmRMSD( tMaskExpr, fit, useMass, debugIn )) return Action::ERR;
  // Initialize reference
  std::string rMaskExpr = actionArgs.GetMaskNext();
  if (rMaskExpr.empty())
    rMaskExpr = tMaskExpr;
  if (REF_.InitRef(previous, first, useMass, fit, reftrajname, REF, RefParm,
                   rMaskExpr, actionArgs, "symmrmsd"))
    return Action::ERR;
  // Set up the RMSD data set. 
  rmsd_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(),"RMSD");
  if (rmsd_==0) return Action::ERR;
  rmsd_->SetScalar( DataSet::M_RMS );
  // Add dataset to data file list
  if (outfile != 0) outfile->AddSet( rmsd_ );
  
  mprintf("    SYMMRMSD: (%s), reference is %s", SRMSD_.TgtMask().MaskString(),
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
Action::RetType Action_SymmetricRmsd::Setup(Topology* currentParm, Topology** parmAddress) {
  // Setup Symmetric RMSD calc (target mask, symmetric atoms etc)
  if (SRMSD_.SetupSymmRMSD( *currentParm )) return Action::ERR;
  // Reference frame setup
  if (REF_.SetupRef(*currentParm, SRMSD_.TgtMask().Nselected(), "symmrmsd"))
    return Action::ERR;
  return Action::OK;
}

// Action_SymmetricRmsd::DoAction()
Action::RetType Action_SymmetricRmsd::DoAction(int frameNum, Frame* currentFrame, 
                                               Frame** frameAddress) 
{
  // Perform any needed reference actions
  REF_.ActionRef( *currentFrame, SRMSD_.Fit(), SRMSD_.UseMass() );
  // Calculate symmetric RMSD
  double rmsdval = SRMSD_.SymmRMSD_CenteredRef( *currentFrame, REF_.RefFrame(), 
                                                REF_.SelectedRef(), REF_.RefTrans() );
  rmsd_->Add(frameNum, &rmsdval);
  if (remap_)
    *frameAddress = (Frame*)SRMSD_.RemapFrame();
  if ( SRMSD_.Fit() )
    (*frameAddress)->Trans_Rot_Trans( SRMSD_.TgtTrans(), SRMSD_.RotMatrix(), REF_.RefTrans() );

  return Action::OK;
}
