// DISTRMSD
#include "Action_DistRmsd.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_DistRmsd::Action_DistRmsd() : drmsd_(0) {}

void Action_DistRmsd::Help() const {
  mprintf("\t[<name>] [<mask>] [<refmask>] [out filename]\n"
          "\t[ first | %s |\n"
          "\t  reftraj <filename> [parm <parmname> | parmindex <#>] ]\n"
          "  Calculate distance RMSD (DME) for specified atoms.\n", DataSetList::RefArgs);
}

// Action_DistRmsd::Init()
/** Called once before traj processing. Set up reference info. */
Action::RetType Action_DistRmsd::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Check for keywords
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  // Reference keywords
  // TODO: Can these just be put in the InitRef call?
  bool first = actionArgs.hasKey("first");
  ReferenceFrame REF = init.DSL().GetReferenceFrame( actionArgs );
  std::string reftrajname = actionArgs.GetStringKey("reftraj");
  Topology* RefParm = init.DSL().GetTopology( actionArgs );
  // Get the RMS mask string for target 
  std::string mask0 = actionArgs.GetMaskNext();
  TgtMask_.SetMaskString(mask0);
  // Get the RMS mask string for reference
  std::string mask1 = actionArgs.GetMaskNext();
  if (mask1.empty())
    mask1 = mask0;

  // Initialize reference
  if (refHolder_.InitRef(false, first, false, false, reftrajname, REF, RefParm,
                         mask1, actionArgs, "distrmsd"))
    return Action::ERR;
 
  // Set up the RMSD data set
  drmsd_ = init.DSL().AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(),"DRMSD");
  if (drmsd_==0) return Action::ERR;
  // Add dataset to data file list
  if (outfile != 0) outfile->AddDataSet( drmsd_ );

  mprintf("    DISTRMSD: (%s), reference is %s\n",TgtMask_.MaskString(),
          refHolder_.RefModeString());

  return Action::OK;
}

// Action_DistRmsd::setup()
/** Called every time the trajectory changes. Set up TgtMask for the new 
  * parmtop and allocate space for selected atoms from the Frame.
  */
Action::RetType Action_DistRmsd::Setup(ActionSetup& setup) {

  if ( setup.Top().SetupIntegerMask(TgtMask_) ) return Action::ERR;
  if ( TgtMask_.None() ) {
    mprintf("Warning: No atoms in mask.\n");
    return Action::SKIP;
  }
  // Allocate space for selected atoms in the frame. This will also put the
  // correct masses in based on the mask.
  SelectedTgt_.SetupFrameFromMask(TgtMask_, setup.Top().Atoms());

  if (refHolder_.SetupRef(setup.Top(), TgtMask_.Nselected(), "distrmsd"))
    return Action::ERR; 

  return Action::OK;
}

// Action_DistRmsd::action()
/** Called every time a frame is read in. Calc distance RMSD.
  * If first is true, set the first frame read in as reference.
  */
Action::RetType Action_DistRmsd::DoAction(int frameNum, ActionFrame& frm) {
  // Perform any needed reference actions
  refHolder_.ActionRef( frm.Frm(), false, false );
  // Set selected frame atoms. Masses have already been set.
  SelectedTgt_.SetCoordinates(frm.Frm(), TgtMask_);
  double DR = SelectedTgt_.DISTRMSD( refHolder_.SelectedRef() );
  drmsd_->Add(frameNum, &DR);
  return Action::OK;
}
