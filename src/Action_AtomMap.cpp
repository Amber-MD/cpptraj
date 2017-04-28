#include <list>
#include "Action_AtomMap.h"
#include "CpptrajStdio.h"
#include "StructureMapper.h"

// CONSTRUCTOR
Action_AtomMap::Action_AtomMap() :
  TgtFrame_(0),
  RefFrame_(0),
  debug_(0),
  newFrame_(0),
  newParm_(0),
  mode_(ALL),
  maponly_(false),
  rmsfit_(false),
  rmsdata_(0)
{}

void Action_AtomMap::Help() const {
  mprintf("\t<target> <reference> [mapout <filename>] [maponly]\n"
          "\t[rmsfit [ rmsout <rmsout> ]] [mode {all | byres}]\n"
          "  Attempt to create a map from atoms in <target> to atoms in <reference> solely\n"
          "  based on how they are bonded; remap <target> so it matches <reference>.\n"
          "  If 'rmsfit' is specified, calculate the RMSD of remapped target to reference.\n"
          "  Modes: 'all'   : try to map all atoms at once (default).\n"
          "         'byres' : map residue by residue (assumes 1 to 1 residue correspondence).\n");
}

// DESTRUCTOR
Action_AtomMap::~Action_AtomMap() {
  if (newFrame_!=0) delete newFrame_;
  if (newParm_!=0) delete newParm_;
}

// Action_AtomMap::Init()
Action::RetType Action_AtomMap::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn; 
  // Get Args
  CpptrajFile* outputfile = init.DFL().AddCpptrajFile(actionArgs.GetStringKey("mapout"),
                                                      "Atom Map",
                                                      DataFileList::TEXT, true);
# ifdef MPI
  // Prevent non-master from writing.
  if (!init.TrajComm().Master()) outputfile = 0;
# endif
  maponly_ = actionArgs.hasKey("maponly");
  rmsfit_ = actionArgs.hasKey("rmsfit");
  std::string modestring = actionArgs.GetStringKey("mode");
  if (!modestring.empty()) {
    if      (modestring == "all") mode_ = ALL;
    else if (modestring == "byres") mode_ = BY_RES;
    else {
      mprinterr("Error: Unrecognized mode: %s\n", modestring.c_str());
      return Action::ERR;
    }
  }
  DataFile* rmsout = 0;
  if (rmsfit_)
    rmsout = init.DFL().AddDataFile( actionArgs.GetStringKey("rmsout"), actionArgs );
  std::string targetName = actionArgs.GetStringNext();
  std::string refName = actionArgs.GetStringNext();
  if (targetName.empty()) {
    mprinterr("Error: No target specified.\n");
    return Action::ERR;
  }
  if (refName.empty()) {
    mprinterr("Error: No reference specified.\n");
    return Action::ERR;
  }
  // Get Reference
  RefFrame_ = (DataSet_Coords_REF*)init.DSL().FindSetOfType( refName, DataSet::REF_FRAME );
  if (RefFrame_ == 0) {
    mprinterr("Error: Could not get reference frame %s\n",refName.c_str());
    return Action::ERR;
  }
  // Get Target
  TgtFrame_ = (DataSet_Coords_REF*)init.DSL().FindSetOfType( targetName, DataSet::REF_FRAME );
  if (TgtFrame_ == 0) {
    mprinterr("Error: Could not get target frame %s\n",targetName.c_str());
    return Action::ERR;
  }
 
  mprintf("    ATOMMAP: Mapping atoms in target topology to given reference.\n"
          "\tTarget topology: '%s'\n\tReference topology: '%s'\n",
          TgtFrame_->Top().c_str(), RefFrame_->Top().c_str());
  if (outputfile != 0)
    mprintf("\tMap will be written to '%s'\n", outputfile->Filename().full());
  if (maponly_)
    mprintf("\tMap will only be written, not used to remap input trajectories.\n");
  else
    mprintf("\tAtoms in input trajectories matching target will be remapped.\n");
  if (!maponly_ && rmsfit_) {
    mprintf("\tWill RMS-fit mapped atoms in tgt to reference.\n");
    if (rmsout != 0) {
      rmsdata_ = init.DSL().AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), "RMSD");
      if (rmsdata_==0) return Action::ERR;
      rmsout->AddDataSet(rmsdata_);
      mprintf("\tRMSDs will be written to '%s'\n", rmsout->DataFilename().full());
    }
  }
  switch (mode_) {
    case BY_RES:
      mprintf("\tCreating map residue-by-residue; assumes 1-to-1 residue correspondence.\n"); break;
    case ALL:
      mprintf("\tCreating map using all atoms at once.\n"); break;
  }

  // ---------------------------------------------
  int err = 0;
  StructureMapper Mapper;
  switch (mode_) {
    case BY_RES: err = Mapper.CreateMapByResidue( RefFrame_, TgtFrame_, debug_ ); break;
    case ALL   : err = Mapper.CreateMap( RefFrame_, TgtFrame_, debug_ ); break;
  }
  if (err != 0) return Action::ERR;
  AMap_ = Mapper.Map();

  // Print atom map
  if (outputfile != 0) {
    Topology const& refTop = RefFrame_->Top();
    Topology const& tgtTop = TgtFrame_->Top();
    outputfile->Printf("%-6s %4s %6s %4s\n","#TgtAt","Tgt","RefAt","Ref");
    for (int refatom = 0; refatom != (int)AMap_.size(); ++refatom) {
      int tgtatom = AMap_[refatom];
      if (tgtatom < 0)
        outputfile->Printf("%6s %4s %6i %4s\n", "---", "---",
                           refatom+1, refTop[refatom].c_str());
      else
        outputfile->Printf("%6i %4s %6i %4s\n",
                           tgtatom+1, tgtTop[tgtatom].c_str(),
                           refatom+1, refTop[refatom].c_str());
    }
  }
  if (maponly_) return Action::OK;

  // If rmsfit is specified, an rms fit of target to reference will be
  // performed using all atoms that were successfully mapped from 
  // target to reference.
  if (rmsfit_) {
    rmsRefFrame_.SetupFrame( Mapper.Nmapped() );
    rmsTgtFrame_ = rmsRefFrame_;
    // Set up a reference frame containing only mapped reference atoms
    rmsRefFrame_.StripUnmappedAtoms(RefFrame_->RefFrame(), AMap_);
    mprintf("\trmsfit: Will rms fit %i atoms from target to reference.\n",
            Mapper.Nmapped());
    return Action::OK;
  }

  // Check if not all atoms could be mapped
  if (Mapper.Nmapped() != (int)AMap_.size()) {
    // If the number of mapped atoms is less than the number of reference
    // atoms but equal to the number of target atoms, can modify the reference
    // frame to only include mapped atoms
    if (Mapper.Nmapped() < (int)AMap_.size() && Mapper.AllTgtMapped()) {
      // Create mask that includes only reference atoms that could be mapped
      AtomMask M1;
      for (int refatom = 0; refatom != (int)AMap_.size(); ++refatom) {
        if (AMap_[refatom] != -1) M1.AddAtom(refatom);
      }
      // Strip reference parm
      mprintf("Warning: Modifying reference '%s' topology and frame to match mapped atoms.\n",
              RefFrame_->legend());
      if (RefFrame_->StripRef( M1 )) return Action::ERR;
      // Since AMap[ ref ] = tgt but ref is now missing any stripped atoms,
      // the indices of AMap must be shifted to match
      int refIndex = 0; // The new index
      for (int refatom = 0; refatom != (int)AMap_.size(); ++refatom) {
        int targetatom = AMap_[refatom];
        if (targetatom > -1)
          AMap_[refIndex++] = targetatom;
      }
    } else {
      mprintf("Warning: Not all atoms were mapped. Frames will not be modified.\n");
      maponly_ = true;
    }
  }

  if (!maponly_) {
    // Set up new Frame
    newFrame_ = new Frame();
    newFrame_->SetupFrameM( TgtFrame_->Top().Atoms() );
    // Set up new Parm
    newParm_ = TgtFrame_->Top().ModifyByMap(AMap_);
  }

  return Action::OK;
}

// Action_AtomMap::Setup()
/** If the current parm does not match the target parm, deactivate. Otherwise
  * replace current parm with mapped parm.
  */
Action::RetType Action_AtomMap::Setup(ActionSetup& setup) {
  if (maponly_) {
    mprintf("\tmaponly was specified, not using atom map during traj read.\n");
    return Action::OK;
  }
  if (setup.Top().Pindex() != TgtFrame_->Top().Pindex() ||
      setup.Top().Natom() != TgtFrame_->Top().Natom()) 
  {
    mprintf("Warning: Map for topology %s -> %s (%i atom).\n",TgtFrame_->Top().c_str(),
            RefFrame_->Top().c_str(), TgtFrame_->Top().Natom());
    mprintf("Warning: Current topology %s (%i atom).\n",setup.Top().c_str(),
            setup.Top().Natom());
    mprintf("Warning: Not using map for this topology.\n");
    return Action::SKIP;
  }
  if (rmsfit_) {
    mprintf("\trmsfit specified, %i atoms.\n",rmsRefFrame_.Natom());
    return Action::OK;
  }
  mprintf("\tMap for parm %s -> %s (%i atom).\n",TgtFrame_->Top().c_str(),
          RefFrame_->Top().c_str(), TgtFrame_->Top().Natom());

  setup.SetTopology( newParm_ );
  
  return Action::MODIFY_TOPOLOGY;
}

// Action_AtomMap::DoAction()
/** Modify the current frame based on the atom map. 
  */
Action::RetType Action_AtomMap::DoAction(int frameNum, ActionFrame& frm) {
  if (maponly_) return Action::OK;

  // Perform RMS fit on mapped atoms only
  if (rmsfit_) {
    // Set target frame up according to atom map.
    rmsTgtFrame_.ModifyByMap(frm.Frm(), AMap_);
    Matrix_3x3 Rot;
    Vec3 Trans, refTrans;
    double R = rmsTgtFrame_.RMSD(rmsRefFrame_, Rot, Trans, refTrans, false);
    frm.ModifyFrm().Trans_Rot_Trans(Trans, Rot, refTrans);
    if (rmsdata_!=0)
      rmsdata_->Add(frameNum, &R);
    return Action::OK;
  }

  // Modify the current frame
  // TODO: Fix this since its probably busted for unmapped atoms
  newFrame_->SetCoordinatesByMap(frm.Frm(), AMap_);
  frm.SetFrame( newFrame_ );
  return Action::MODIFY_COORDS;
}
