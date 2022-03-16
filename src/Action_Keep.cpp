#include "Action_Keep.h"
#include "CharMask.h"
#include "CpptrajStdio.h"
#include "DataSet_string.h"
#include "Range.h"

/** CONSTRUCTOR */
Action_Keep::Action_Keep() :
  debug_(0),
  currentParm_(0),
  keepParm_(0),
  bridgeData_(0),
  nbridge_(0),
  bridgeWarn_(true),
  nNonBridgeAtoms_(-1)
{}

/** DESTRUCTOR */
Action_Keep::~Action_Keep() {
  if (keepParm_ != 0) delete keepParm_;
}

// Residue status characters
const char Action_Keep::STAT_NONE_ = 'X';
const char Action_Keep::STAT_BRIDGERES_ = 'B';
const char Action_Keep::STAT_NONBRIDGERES_ = 'U';

// Action_Keep::Help()
void Action_Keep::Help() const {
  mprintf("\t[ bridgedata <bridge data set> [nbridge <#>] [nobridgewarn]\n"
          "\t [bridgeresname <res name>] bridgeresonly <resrange>] ]\n"
          "\t[keepmask <atoms to keep>]\n");
  mprintf("%s", ActionTopWriter::Keywords());
  mprintf("  Keep only specified parts of the system.\n");
  mprintf("%s", ActionTopWriter::Options());
}

// Action_Keep::Init()
Action::RetType Action_Keep::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  std::string keepmaskstr = actionArgs.GetStringKey("keepmask");
  if (!keepmaskstr.empty()) {
    if ( keepMask_.SetMaskString( keepmaskstr )) {
      mprinterr("Error: Invalid mask for 'keepmask'\n");
      return Action::ERR;
    }
  }

  std::string bridgeDataName = actionArgs.GetStringKey("bridgedata");
  if (!bridgeDataName.empty()) {
    DataSet* ds = init.DSL().GetDataSet( bridgeDataName );
    if (ds == 0) {
      mprinterr("Error: No data set corresponding to '%s'\n", bridgeDataName.c_str());
      return Action::ERR;
    }
    if (ds->Type() != DataSet::STRING) {
      mprinterr("Error: Bridge data set '%s' is not a string data set.\n", ds->legend());
      return Action::ERR;
    }
    bridgeData_ = static_cast<DataSet_string*>( ds );
    nbridge_ = actionArgs.getKeyInt("nbridge", 1);
    if (nbridge_ < 1) {
      mprinterr("Error: Number of bridging residues to keep must be >= 1.\n");
      return Action::ERR;
    }
    bridgeResName_ = actionArgs.GetStringKey("bridgeresname", "WAT");
    std::string bridgeresonlystr = actionArgs.GetStringKey("bridgeresonly");
    if (!bridgeresonlystr.empty()) {
      Range brange;
      if (brange.SetRange( bridgeresonlystr )) {
        mprinterr("Error: Invalid range given for 'bridgeresonly': %s\n", bridgeresonlystr.c_str());
        return Action::ERR;
      }
      // User residue numbers start from 1
      bridgeResOnly_.clear();
      for (Range::const_iterator it = brange.begin(); it != brange.end(); ++it)
        bridgeResOnly_.push_back( *it - 1 );
    }
    bridgeWarn_ = !actionArgs.hasKey("nobridgewarn");
  } else if (!keepMask_.MaskStringSet()) {
    mprinterr("Error: Nothing specified to keep.\n");
    return Action::ERR;
  }

  topWriter_.InitTopWriter(actionArgs, "keep", debugIn);

  mprintf("    KEEP:\n");
  if (keepMask_.MaskStringSet()) {
    mprintf("\tKeeping atoms selected by '%s'\n", keepMask_.MaskString());
  }
  if (bridgeData_ != 0) {
    mprintf("\tKeeping bridging residues from bridge ID data set: %s\n", bridgeData_->legend());
    mprintf("\t# of bridging residues to keep: %i\n", nbridge_);
    mprintf("\tBridge residue name: %s\n", bridgeResName_.c_str());
    if (!bridgeResOnly_.empty()) {
      mprintf("\tOnly keeping bridge residues when bridging residues:");
      for (Iarray::const_iterator it = bridgeResOnly_.begin(); it != bridgeResOnly_.end(); ++it)
        mprintf(" %i", *it + 1);
      mprintf("\n");
    }
    if (bridgeWarn_)
      mprintf("\tWill warn when # active bridges does not match requested.\n");
    else
      mprintf("\tHiding warnings for when # active bridges does not match requested.\n");
#   ifdef MPI
    if (init.TrajComm().Size() > 1)
      mprintf("Warning: Skipping frames with incorrect # of bridges in parallel\n"
              "Warning:   can cause certain actions (e.g. 'rms') to hang.\n"
              "Warning:   In addition, trajectories written after skipping frames may have issues.\n");
#   endif
  }
  topWriter_.PrintOptions();
  return Action::OK;
}

// Action_Keep::Setup()
Action::RetType Action_Keep::Setup(ActionSetup& setup)
{
  currentParm_ = setup.TopAddress();

  atomsToKeep_.ClearSelected();
  atomsToKeep_.SetNatoms( setup.Top().Natom() );

  resStat_.assign( setup.Top().Nres(), STAT_NONE_ );

  if (keepMask_.MaskStringSet()) {
    if (setup.Top().SetupIntegerMask( keepMask_ )) {
      mprinterr("Error: Could not set up keep mask '%s'\n", keepMask_.MaskString());
      return Action::ERR;
    }
    if (keepMask_.None()) {
      mprintf("Warning: No atoms selected for keep mask.\n");
      return Action::SKIP;
    }
    keepMask_.MaskInfo();
    for (AtomMask::const_iterator it = keepMask_.begin(); it != keepMask_.end(); ++it)
      atomsToKeep_.AddSelectedAtom(*it);
  }

  if (bridgeData_ != 0) {
    // Set up to keep bridge residues
    AtomMask bmask;
    if (bmask.SetMaskString( ":" + bridgeResName_ )) {
      mprinterr("Error: Could not set up mask for bridge residues: %s\n", bridgeResName_.c_str());
      return Action::ERR;
    }
    // Select potential bridge residues
    if (setup.Top().SetupIntegerMask( bmask )) {
      mprinterr("Error: Setting up bridge residue mask failed.\n");
      return Action::ERR;
    }
    if (bmask.None()) {
      mprintf("Warning: No potential bridge residues selected.\n");
      return Action::SKIP;
    }
    bmask.MaskInfo();
    // Ensure all bridge residues have same # atoms
    std::vector<int> Rnums = setup.Top().ResnumsSelectedBy( bmask );
    int resSize = -1;
    for (std::vector<int>::const_iterator rnum = Rnums.begin();
                                          rnum != Rnums.end(); ++rnum)
    {
      if (resSize == -1)
        resSize = setup.Top().Res(*rnum).NumAtoms();
      else if (setup.Top().Res(*rnum).NumAtoms() != resSize) {
        mprinterr("Error: Residue '%s' size (%i) != first residue size (%s)\n",
                  setup.Top().TruncResNameNum(*rnum).c_str(),
                  setup.Top().Res(*rnum).NumAtoms(), resSize);
        return Action::ERR;
      }
      resStat_[*rnum] = STAT_BRIDGERES_;
    }
    mprintf("\tBridge residue size= %i\n", resSize);
    // This character mask will be used if keepMask is not set.
    CharMask cmask;
    if (!keepMask_.MaskStringSet()) {
      // Keep all atoms not in the bridge mask
      cmask = CharMask( bmask.ConvertToCharMask(), bmask.Nselected() );
      for (int idx = 0; idx != setup.Top().Natom(); idx++) {
        if (!cmask.AtomInCharMask(idx)) {
          atomsToKeep_.AddSelectedAtom( idx );
        }
      }
    }
    // Set up non-bridge residue status
    if (!bridgeResOnly_.empty()) {
      for (Iarray::const_iterator it = bridgeResOnly_.begin(); it != bridgeResOnly_.end(); ++it)
      {
        if (*it < 0 || *it >= setup.Top().Nres()) {
          mprinterr("Error: Residue # %i for bridgeresonly out of range.\n", *it+1);
          return Action::ERR;
        }
        resStat_[ *it ] = STAT_NONBRIDGERES_;
      }
    } else if (keepMask_.MaskStringSet()) {
      for (AtomMask::const_iterator at = keepMask_.begin(); at != keepMask_.end(); ++at)
        resStat_[ setup.Top()[*at].ResNum() ] = STAT_NONBRIDGERES_;
    } else {
      for (int idx = 0; idx != setup.Top().Natom(); idx++) {
        if (!cmask.AtomInCharMask(idx))
          resStat_[ setup.Top()[idx].ResNum() ] = STAT_NONBRIDGERES_;
      }
    }
    // Will keep only the first nbridge_ residues
    if ((unsigned int)nbridge_ > Rnums.size()) {
      mprinterr("Error: Number of bridge residues to keep (%i) > number of potential bridge residues (%zu).\n", nbridge_, Rnums.size());
      return Action::ERR;
    }
    nNonBridgeAtoms_ = atomsToKeep_.Nselected();
    for (int ridx = 0; ridx < nbridge_; ridx++) {
      Residue const& res = setup.Top().Res(Rnums[ridx]);
      for (int at = res.FirstAtom(); at != res.LastAtom(); at++) {
        atomsToKeep_.AddSelectedAtom( at );
      }
    }
  } // END bridgeData

  // DEBUG: Print res stat
  if (debug_ > 1) {
    for (int rnum = 0; rnum != setup.Top().Nres(); rnum++) {
      mprintf("DEBUG: Res %20s stat %c\n", setup.Top().TruncResNameNum(rnum).c_str(), resStat_[rnum]);
    }
  }

  // Create topology with only atoms to keep
  if (keepParm_ != 0)
    delete keepParm_;
  if (atomsToKeep_.Nselected() < 1) {
    mprintf("Warning: No atoms to keep.\n");
    return Action::SKIP;
  }
  keepParm_ = setup.Top().modifyStateByMask( atomsToKeep_ );
  if (keepParm_ == 0) {
    mprinterr("Error: Could not create topology for kept atoms.\n");
    return Action::ERR;
  }
  setup.SetTopology( keepParm_ );
  // Modify State if asked
  if (topWriter_.ModifyActionState(setup, keepParm_))
    return Action::ERR;
  keepParm_->Brief("Topology for kept atoms:");
  keepFrame_.SetupFrameV( setup.Top().Atoms(), setup.CoordInfo() );

  topWriter_.WriteTops( *keepParm_ );
    
  return Action::MODIFY_TOPOLOGY;
}

// Action_Keep::DoAction()
Action::RetType Action_Keep::DoAction(int frameNum, ActionFrame& frm)
{
  Action::RetType err = Action::OK;
  if (bridgeData_ != 0)
    err = keepBridge(frameNum, frm);

  if (err == Action::OK) {
    keepFrame_.SetFrame(frm.Frm(), atomsToKeep_);
    frm.SetFrame( &keepFrame_ );
  }

  return err;
}

/** Want to keep only residues specified in a bridge ID data set. */
Action::RetType Action_Keep::keepBridge(int frameNum, ActionFrame& frm) {
  atomsToKeep_.ShrinkSelectedTo( nNonBridgeAtoms_ );
  // Ensure we can get data
  if ((unsigned int)frameNum >= bridgeData_->Size()) {
    mprinterr("Error: Frame # %i is out of range for bridge data '%s' (size is %zu)\n",
              frameNum+1, bridgeData_->legend(), bridgeData_->Size());
    return Action::ERR;
  }
  std::string const& bridgeIDstr = (*bridgeData_)[frameNum];
  if (bridgeIDstr == "None") {
    if (bridgeWarn_)
      mprintf("Warning: Frame %i has no bridging residues.\n", frameNum+1);
    return Action::SUPPRESS_COORD_OUTPUT;
  }
  ArgList bridgeID( bridgeIDstr, "," );
  if (debug_ > 0)
    mprintf("DEBUG: Frame %i has %i bridging residues.\n", frameNum+1, bridgeID.Nargs());
  if (bridgeID.Nargs() < nbridge_) {
    if (bridgeWarn_)
      mprintf("Warning: Frame %i has fewer total bridges than requested (%i).\n", frameNum+1, bridgeID.Nargs());
    return Action::SUPPRESS_COORD_OUTPUT;
  }
  int nBridgeInFrame = 0;
  for (int nb = 0; nb != bridgeID.Nargs(); nb++) {
    // Format: <bres#>(ures0+ures1+...)
    ArgList bridge( bridgeID[nb], "()+" );
    if (bridge.Nargs() < 3) {
      mprinterr("Error: Expected at least 3 args for bridge ID '%s', got %i\n",
                bridgeID[nb].c_str(), bridge.Nargs());
      return Action::ERR;
    }
    // bres# indices start from 1
    int bres = bridge.getNextInteger(0) - 1;
    if (debug_ > 0)
      mprintf("DEBUG:   Bridge res %i\n", bres+1);
    // Check that bres is actually a bridging residue
    if (bres < 0) {
      mprinterr("Error: Invalid bridging residue # %i for bridge '%s'\n", bres+1, bridge.ArgLine());
    } else if ( resStat_[bres] != STAT_BRIDGERES_ ) {
      mprinterr("Error: Residue %s listed as bridging but was not selected by 'bridgeresname'.\n",
                currentParm_->TruncResNameNum(bres).c_str());
    }
    // Check that residues being bridged are active
    bool bridgeIsActive = true;
    int nBres = bridge.getNextInteger(0) - 1;
    while (nBres > -1) {
      if ( resStat_[nBres] != STAT_NONBRIDGERES_ ) {
        bridgeIsActive = false;
        break;
      }
      nBres = bridge.getNextInteger(0) - 1;
    }
    if (bridgeIsActive) {
      nBridgeInFrame++;
      if (nBridgeInFrame > nbridge_) {
        if (bridgeWarn_)
          mprintf("Warning: More active bridges in frame %i (%i) than specified (%i); skipping.\n",
                  frameNum+1, nBridgeInFrame, nbridge_);
        return SUPPRESS_COORD_OUTPUT;
      }
      Residue const& res = currentParm_->Res(bres);
      for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
        atomsToKeep_.AddSelectedAtom( at );
    }
  }
  if (nBridgeInFrame < nbridge_) {
    if (bridgeWarn_)
      mprintf("Warning: Fewer active bridges in frame %i (%i) than specified (%i); skipping.\n",
              frameNum+1, nBridgeInFrame, nbridge_);
    return SUPPRESS_COORD_OUTPUT;
  }

  return Action::OK;
}
