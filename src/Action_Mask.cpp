#include "Action_Mask.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Mask::Action_Mask() :
  outfile_(0),
  fnum_(0),
  anum_(0),
  aname_(0),
  rnum_(0),
  rname_(0),
  mnum_(0),
  idx_(0),
  CurrentParm_(0),
  debug_(0),
  writeTraj_(false)
{} 

void Action_Mask::Help() const {
  mprintf("\t<mask1> [maskout <filename>] [maskpdb <filename> | maskmol2 <filename>]\n"
          "  Print atoms selected by <mask1> to file specified by 'maskout' and/or\n"
          "  the PDB or Mol2 file specified by 'maskpdb' or 'maskmol2'. Good for\n"
          "  distance-based masks.\n");
}

// Action_Mask::Init()
// NOTE: Could also split the arglist at maskpdb and make it so any type of 
//       file can be written out.
Action::RetType Action_Mask::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  // Get Keywords
  outfile_ = init.DFL().AddCpptrajFile(actionArgs.GetStringKey("maskout"), "Atoms in mask");
# ifdef MPI
  trajComm_ = init.TrajComm();
  if ( trajComm_.Size() > 1 && outfile_ != 0) {
    mprinterr("Error: 'maskout' currently only works with 1 process (currently %i)\n"
              "Error:   Consider using the 'name'/'out' keywords instead\n", trajComm_.Size());
    return Action::ERR;
  }
# endif
  std::string maskpdb = actionArgs.GetStringKey("maskpdb");
  std::string maskmol2 = actionArgs.GetStringKey("maskmol2");
  std::string dsname = actionArgs.GetStringKey("name");
  std::string dsout = actionArgs.GetStringKey("out");
  // At least 1 of maskout, maskpdb, maskmol2, or name must be specified.
  if (outfile_ == 0 && maskpdb.empty() && maskmol2.empty() && dsname.empty()) {
    mprinterr("Error: At least one of maskout, maskpdb, maskmol2, or name must be specified.\n");
    return Action::ERR;
  }
  // Set up any trajectory options
  TrajectoryFile::TrajFormatType trajFmt = TrajectoryFile::PDBFILE;
  if (!maskpdb.empty() || !maskmol2.empty()) {
    outtraj_.SetDebug( debug_ );
    ArgList trajArgs;
    if (!maskpdb.empty()) {
      // Set pdb output options: multi so that 1 file per frame is written; dumpq
      // so that charges are written out.
      trajArgs = ArgList("multi dumpq nobox");
    } else if (!maskmol2.empty()) {
      maskpdb = maskmol2;
      trajFmt = TrajectoryFile::MOL2FILE;
      trajArgs = ArgList("multi nobox");
    }
    if (outtraj_.InitEnsembleTrajWrite(maskpdb, trajArgs, trajFmt, init.DSL().EnsembleNum()))
      return Action::ERR;
    writeTraj_ = true;
  } else
    writeTraj_ = false;
  // Get Mask
  if (Mask1_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;
  // Set up data sets
  if (!dsname.empty() || !dsout.empty()) {
    if (dsname.empty()) dsname = init.DSL().GenerateDefaultName("MASK");
    MetaData::tsType ts = MetaData::NOT_TS; // None are a straight time series
    fnum_  = init.DSL().AddSet(DataSet::INTEGER, MetaData(dsname, "Frm",   ts));
    anum_  = init.DSL().AddSet(DataSet::INTEGER, MetaData(dsname, "AtNum", ts));
    aname_ = init.DSL().AddSet(DataSet::STRING,  MetaData(dsname, "Aname", ts));
    rnum_  = init.DSL().AddSet(DataSet::INTEGER, MetaData(dsname, "Rnum",  ts));
    rname_ = init.DSL().AddSet(DataSet::STRING,  MetaData(dsname, "Rname", ts));
    mnum_  = init.DSL().AddSet(DataSet::INTEGER, MetaData(dsname, "Mnum",  ts));
    if (fnum_ == 0 || anum_ == 0 || aname_ == 0 ||
        rnum_ == 0 || rname_ == 0 || mnum_ == 0)
      return Action::ERR;
    DataFile* dout = init.DFL().AddDataFile( dsout, "noxcol", actionArgs );
    if (dout != 0) {
      dout->AddDataSet( fnum_ );
      dout->AddDataSet( anum_ );
      dout->AddDataSet( aname_ );
      dout->AddDataSet( rnum_ );
      dout->AddDataSet( rname_ );
      dout->AddDataSet( mnum_ );
    }
    idx_ = 0;
  }
  mprintf("    ACTIONMASK: Information on atoms in mask %s will be printed",
          Mask1_.MaskString());
  if (outfile_ != 0)
    mprintf(" to file %s",outfile_->Filename().full());
  mprintf(".\n");
  if (writeTraj_) 
    mprintf("\t%ss of atoms in mask will be written to %s.X\n",
            TrajectoryFile::FormatString(trajFmt), outtraj_.Traj().Filename().full());
  if (fnum_ != 0)
    mprintf("\tData sets will be saved with name '%s'\n", fnum_->Meta().Name().c_str());
  // Header
  if (outfile_ != 0)
    outfile_->Printf("%-8s %8s %4s %8s %4s %8s\n","#Frame","AtomNum","Atom",
                     "ResNum","Res", "MolNum");
  return Action::OK;
}

// Action_Mask::Setup()
Action::RetType Action_Mask::Setup(ActionSetup& setup) {
  CurrentParm_ = setup.TopAddress();
  cInfo_ = setup.CoordInfo();
  return Action::OK;
}

// Action_Mask::DoAction()
Action::RetType Action_Mask::DoAction(int frameNum, ActionFrame& frm) {
  // Get atom selection
  if ( CurrentParm_->SetupCharMask(Mask1_, frm.Frm()) ) {
    mprintf("Warning: Could not set up atom mask [%s]\n",
            Mask1_.MaskString());
    return Action::ERR;
  }
  // Print out information for every atom in the mask
  for (int atom=0; atom < CurrentParm_->Natom(); atom++) {
    if (Mask1_.AtomInCharMask(atom)) {
      int res = (*CurrentParm_)[atom].ResNum();
      int fn = frm.TrajoutNum() + 1;
      int an = atom + 1;
      int rn = res + 1;
      int mn = (*CurrentParm_)[atom].MolNum() + 1;
      if (outfile_ != 0)
        outfile_->Printf("%8i %8i %4s %8i %4s %8i\n", fn, an, (*CurrentParm_)[atom].c_str(),
                         rn, CurrentParm_->Res(res).c_str(), mn);
      /*mprintf(" Type=%4s",CurrentParm_->types[atom]);
      mprintf(" Charge=%lf",CurrentParm_->charge[atom]);
      mprintf(" Mass=%lf",CurrentParm_->mass[atom]);
      outfile.Printf("\n");*/
      if (fnum_ != 0) {
        fnum_->Add( idx_, &fn );
        anum_->Add( idx_, &an );
        rnum_->Add( idx_, &rn );
        mnum_->Add( idx_, &mn );
        aname_->Add( idx_, (*CurrentParm_)[atom].c_str() );
        rname_->Add( idx_, CurrentParm_->Res(res).c_str() );
        idx_++;
      }
    }
  }

  // Optional write out of selected atoms for the frame.
  if (writeTraj_) {
    // Convert Mask1 to an integer mask for use in parm/frame functions
    AtomMask Mask2( Mask1_.ConvertToIntMask(), Mask1_.Natom() );
    // Create new parm and frame based on atoms in Mask. Since we dont care
    // about advanced parm info for PDB write just do a partial modify.
    Topology* pdbParm = CurrentParm_->partialModifyStateByMask(Mask2);
    //pdbParm->Summary(); // DEBUG
    Frame pdbFrame(frm.Frm(), Mask2);
    // Set up output trajectory file. 
    if (outtraj_.SetupTrajWrite(pdbParm, cInfo_, 1))
    {
      mprinterr("Error: %s: Could not write mask atoms for frame %i.\n",
                outtraj_.Traj().Filename().base(), frm.TrajoutNum() + 1);
    } else {
      if (debug_ > 0) outtraj_.PrintInfo(0);
      outtraj_.WriteSingle(frm.TrajoutNum(), pdbFrame);
      outtraj_.EndTraj();
    }
    delete pdbParm;
  }

  return Action::OK;
}

#ifdef MPI
/** Since datasets are actually # frames * however many atoms found in mask
  * at each frame, sync here.
  */
int Action_Mask::SyncAction() {
  if (fnum_ == 0) return 0;
  // Get total number of mask entries.
  std::vector<int> rank_frames( trajComm_.Size() );
  trajComm_.GatherMaster( &idx_, 1, MPI_INT, &(rank_frames[0]) );
  for (int rank = 1; rank < trajComm_.Size(); rank++)
    idx_ += rank_frames[ rank ];
  fnum_->Sync( idx_, rank_frames, trajComm_ );
  fnum_->SetNeedsSync( false );
  anum_->Sync( idx_, rank_frames, trajComm_ );
  anum_->SetNeedsSync( false );
  aname_->Sync( idx_, rank_frames, trajComm_ );
  aname_->SetNeedsSync( false );
  rnum_->Sync( idx_, rank_frames, trajComm_ );
  rnum_->SetNeedsSync( false );
  rname_->Sync( idx_, rank_frames, trajComm_ );
  rname_->SetNeedsSync( false );
  mnum_->Sync( idx_, rank_frames, trajComm_ );
  mnum_->SetNeedsSync( false );

  return 0;
}
#endif
