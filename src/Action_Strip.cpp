// Action_Strip
#include "Action_Strip.h"
#include "CpptrajStdio.h"
#include "ParmFile.h"

// CONSTRUCTOR
Action_Strip::Action_Strip() :
  newParm_(0), newCinfo_(0), masterDSL_(0), removeBoxInfo_(false) {}

void Action_Strip::Help() const {
  mprintf("\t<mask> [outprefix <name>] [parmout <file>] [nobox]\n"
          "  Strip atoms in <mask> from the system.\n");
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
  prefix_ = actionArgs.GetStringKey("outprefix");
  parmoutName_ = actionArgs.GetStringKey("parmout");
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
  if (!prefix_.empty()) 
    mprintf("\tStripped topology will be output with prefix '%s'\n", prefix_.c_str());
  if (!parmoutName_.empty())
    mprintf("\tStripped topology will be output with name '%s'\n", parmoutName_.c_str());
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
    mprintf("Warning: strip: Mask [%s] has no atoms.\n",M1_.MaskString());
    return Action::SKIP;
  }
  int numStripped = setup.Top().Natom() - M1_.Nselected();
  mprintf("\tStripping %i atoms.\n", numStripped);
  // If no atoms will be stripped, no need to use this command. SKIP
  if ( numStripped == 0 ) {
    mprintf("Warning: No atoms to strip. Skipping 'strip' for topology '%s'\n",
            setup.Top().c_str());
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
  if (!prefix_.empty()) {
    ParmFile pfile;
    if ( pfile.WritePrefixTopology( setup.Top(), prefix_, ParmFile::AMBERPARM, 0 ) )
      mprinterr("Error: Could not write out stripped topology file.\n");
  }
  if (!parmoutName_.empty()) {
    ParmFile pfile;
    if ( pfile.WriteTopology( setup.Top(), parmoutName_, ParmFile::AMBERPARM, 0 ) )
      mprinterr("Error: Could not write out stripped topology file '%s'\n", parmoutName_.c_str());
  }

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
