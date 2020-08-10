// Action_Strip
#include "Action_Strip.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Strip::Action_Strip() :
  newParm_(0), newCinfo_(0), masterDSL_(0), removeBoxInfo_(false) {}

void Action_Strip::Help() const {
  mprintf("\t<mask> [nobox]\n");
  mprintf("%s", ActionTopWriter::Keywords());
  mprintf("  Strip atoms in <mask> from the system. If 'nobox' specified,\n"
          "  remove any unit cell info.\n");
  mprintf("%s", ActionTopWriter::Options());
}

// DESTRUCTOR
Action_Strip::~Action_Strip() {
  if (newParm_!=0) delete newParm_;
  if (newCinfo_ != 0) delete newCinfo_;
}

// Action_Strip::Init()
Action::RetType Action_Strip::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get output stripped parm filename
  topWriter_.InitTopWriter(actionArgs, "stripped", debugIn);
  removeBoxInfo_ = actionArgs.hasKey("nobox");

  // Get mask of atoms to be stripped
  std::string mask1 = actionArgs.GetMaskNext();
  //mprintf("    Mask 1: %s\n",mask1);
  if (mask1.empty()) {
    mprinterr("Error: Requires atom mask.\n");
    return Action::ERR;
  }
  if (M1_.SetMaskString(mask1)) return Action::ERR;
  // We want to strip the atoms inside the mask and keep those outside
  // the mask. Since modifyStateByMask needs to know the kept atoms,
  // invert the mask selection.
  M1_.InvertMaskExpression();

  mprintf("    STRIP: Stripping atoms in mask [%s]\n", M1_.MaskString());
  topWriter_.PrintOptions();
  if (removeBoxInfo_)
    mprintf("\tAny existing box information will be removed.\n");
  masterDSL_ = init.DslPtr();
  return Action::OK;
}

// Action_Strip::Setup()
/** Attempt to create a new stripped down version of the input parmtop
  */
Action::RetType Action_Strip::Setup(ActionSetup& setup) {
  if (setup.Top().SetupIntegerMask( M1_ )) return Action::ERR;
  if (M1_.None()) {
    // If no atoms will be kept, no need for this command. SKIP.
    mprintf("Warning: Mask [%s] would strip all atoms. Skipping.\n", M1_.MaskString());
    return Action::SKIP;
  }
  int numStripped = setup.Top().Natom() - M1_.Nselected();
  mprintf("\tStripping %i atoms.\n", numStripped);
  // SANITY CHECK: If no atoms will be stripped, no need to use this command. SKIP
  if ( numStripped == 0 ) {
    mprintf("Warning: No atoms to strip. Skipping.\n");
    return Action::SKIP;
  }

  // Attempt to create new parmtop based on mask, set as new Topology
  if (newParm_ != 0) delete newParm_;
  newParm_ = setup.Top().modifyStateByMask(M1_);
  if (newParm_ == 0) {
    mprinterr("Error: Could not create new topology.\n");
    return Action::ERR;
  }
  setup.SetTopology( newParm_ );
  // Remove box information if asked
  if (removeBoxInfo_) {
    newParm_->SetParmBox( Box() );
    newCinfo_ = new CoordinateInfo( setup.CoordInfo() );
    newCinfo_->SetBox( Box() );
    setup.SetCoordInfo( newCinfo_ );
  }
  newParm_->Brief("Stripped topology:");
  // Allocate space for new frame
  newFrame_.SetupFrameV(setup.Top().Atoms(), setup.CoordInfo());

  // If prefix given then output stripped topology
  topWriter_.WriteTops( setup.Top() );

  return Action::MODIFY_TOPOLOGY;
}

// Action_Strip::DoAction()
/** Modify the coordinate frame to reflect stripped parmtop. */
Action::RetType Action_Strip::DoAction(int frameNum, ActionFrame& frm) {

  newFrame_.SetFrame(frm.Frm(), M1_);

  // Set frame
  frm.SetFrame( &newFrame_ );

  return Action::MODIFY_COORDS;
}
