#include "Action_Keep.h"
#include "CharMask.h"
#include "CpptrajStdio.h"
#include "DataSet_string.h"

/** CONSTRUCTOR */
Action_Keep::Action_Keep() :
  currentParm_(0),
  keepParm_(0),
  bridgeData_(0),
  nbridge_(0)
{}

/** DESTRUCTOR */
Action_Keep::~Action_Keep() {
  if (keepParm_ != 0) delete keepParm_;
}

// Action_Keep::Help()
void Action_Keep::Help() const {
  mprintf("\t[bridgedata <bridge data set> [nbridge <#>] [bridgeresname <res name>]\n"
          "\t[keepmask <atoms to keep>]\n");
  mprintf("%s", ActionTopWriter::Keywords());
  mprintf("  Keep only specified parts of the system.\n");
  mprintf("%s", ActionTopWriter::Options());
}

// Action_Keep::Init()
Action::RetType Action_Keep::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
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
  } else {
    mprinterr("Error: Nothing specified to keep.\n");
    return Action::ERR;
  }

  topWriter_.InitTopWriter(actionArgs, "keep", debugIn);

  mprintf("    KEEP:\n");
  if (bridgeData_ != 0) {
    mprintf("\tBridge ID data set: %s\n", bridgeData_->legend());
    mprintf("\t# of bridging residues to keep: %i\n", nbridge_);
    mprintf("\tBridge residue name: %s\n", bridgeResName_.c_str());
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
    }
    mprintf("\tBridge residue size= %i\n", resSize);
    if (!keepMask_.MaskStringSet()) {
      // Keep all atoms not in the bridge mask
      CharMask cmask( bmask.ConvertToCharMask(), bmask.Nselected() );
      for (int idx = 0; idx != setup.Top().Natom(); idx++)
        if (!cmask.AtomInCharMask(idx))
          atomsToKeep_.AddSelectedAtom( idx );
    }
    // Will keep only the first nbridge_ residues
    if ((unsigned int)nbridge_ > Rnums.size()) {
      mprinterr("Error: Number of bridge residues to keep (%i) > number of potential bridge residues (%zu).\n", nbridge_, Rnums.size());
      return Action::ERR;
    }
    //idxMaskPair_.clear();
    nNonBridgeAtoms_ = atomsToKeep_.Nselected();
    for (int ridx = 0; ridx < nbridge_; ridx++) {
      Residue const& res = setup.Top().Res(Rnums[ridx]);
      //idxMaskPair_.push_back( IdxMaskPairType(selectedIdx, AtomMask()) );
      for (int at = res.FirstAtom(); at != res.LastAtom(); at++) {
        atomsToKeep_.AddSelectedAtom( at );
        //idxMaskPair_.back().second.AddSelectedAtom(
      }
    }
    

  } // END bridgeData

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
    mprintf("DEBUG: Frame %i has no bridging residues.\n", frameNum+1);
    return Action::SUPPRESS_COORD_OUTPUT;
  }
  ArgList bridgeID( bridgeIDstr, "," );
  mprintf("DEBUG: Frame %i has %i bridging residues.\n", frameNum+1, bridgeID.Nargs());
  if (bridgeID.Nargs() < nbridge_) {
    mprintf("Warning: Frame %i has fewer bridges than requested (%i).\n", frameNum+1, bridgeID.Nargs());
    return Action::SUPPRESS_COORD_OUTPUT;
  }
  for (int nb = 0; nb != bridgeID.Nargs(); nb++) {
    // Format: <bres#>(ures0+ures1+...)
    // bres# indices start from 1
    ArgList bridge( bridgeID[nb], "()+" );
    if (bridge.Nargs() < 3) {
      mprinterr("Error: Expected at least 3 args for bridge ID '%s', got %i\n",
                bridgeID[nb].c_str(), bridge.Nargs());
      return Action::ERR;
    }
    int bres = bridge.getNextInteger(-1);
    mprintf("DEBUG: Bridge res %i\n", bres);
    // TODO check that bres is in resnums

    Residue const& res = currentParm_->Res(bres-1);
    for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
      atomsToKeep_.AddSelectedAtom( at );
  }

  return Action::OK;
}
