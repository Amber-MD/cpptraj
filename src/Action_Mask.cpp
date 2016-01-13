#include "Action_Mask.h"
#include "CpptrajStdio.h"
#include "Trajout_Single.h"

// CONSTRUCTOR
Action_Mask::Action_Mask() :
  ensembleNum_(-1),
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
  trajFmt_(TrajectoryFile::PDBFILE),
  trajOpt_(0)
{ } 

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
  ensembleNum_ = init.DSL().EnsembleNum();
  debug_ = debugIn;
  // Get Keywords
  outfile_ = init.DFL().AddCpptrajFile(actionArgs.GetStringKey("maskout"), "Atoms in mask");
  maskpdb_ = actionArgs.GetStringKey("maskpdb");
  std::string maskmol2 = actionArgs.GetStringKey("maskmol2");
  std::string dsname = actionArgs.GetStringKey("name");
  std::string dsout = actionArgs.GetStringKey("out");
  // At least 1 of maskout, maskpdb, maskmol2, or name must be specified.
  if (outfile_ == 0 && maskpdb_.empty() && maskmol2.empty() && dsname.empty()) {
    mprinterr("Error: At least one of maskout, maskpdb, maskmol2, or name must be specified.\n");
    return Action::ERR;
  }
  // Set up any trajectory options
  if (!maskpdb_.empty()) {
    trajFmt_ = TrajectoryFile::PDBFILE;
    // Set pdb output options: multi so that 1 file per frame is written; dumpq
    // so that charges are written out.
    trajOpt_ = "multi dumpq nobox";
  } else if (!maskmol2.empty()) {
    maskpdb_ = maskmol2;
    trajFmt_ = TrajectoryFile::MOL2FILE;
    trajOpt_ = "multi nobox";
  }
  // Get Mask
  Mask1_.SetMaskString( actionArgs.GetMaskNext() );
  // Set up data sets
  if (!dsname.empty()) {
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
    DataFile* dout = init.DFL().AddDataFile( dsout, actionArgs );
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
  if (!maskpdb_.empty()) 
    mprintf("\t%ss of atoms in mask will be written to %s.X\n",
            TrajectoryFile::FormatString(trajFmt_), maskpdb_.c_str());
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
  currentCoordInfo_ = setup.CoordInfo();
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
  if (!maskpdb_.empty()) {
    Trajout_Single coordsOut;
    // Convert Mask1 to an integer mask for use in parm/frame functions
    AtomMask Mask2( Mask1_.ConvertToIntMask(), Mask1_.Natom() );
    // Create new parm and frame based on atoms in Mask. Since we dont care
    // about advanced parm info for PDB write just do a partial modify.
    Topology* pdbParm = CurrentParm_->partialModifyStateByMask(Mask2);
    //pdbParm->Summary(); // DEBUG
    Frame pdbFrame(frm.Frm(), Mask2);
    // Set up output trajectory file. 
    coordsOut.SetDebug(debug_);
    if (coordsOut.PrepareEnsembleTrajWrite(maskpdb_,trajOpt_,pdbParm,
                                           currentCoordInfo_,
                                           1,trajFmt_,ensembleNum_)) 
    {
      mprinterr("Error: %s: Could not write mask atoms for frame %i.\n",
                maskpdb_.c_str(), frameNum + 1);
    } else {
      if (debug_ > 0) coordsOut.PrintInfo(0);
      coordsOut.WriteSingle(frameNum, pdbFrame);
      coordsOut.EndTraj();
    }
    delete pdbParm;
  }

  return Action::OK;
}

#ifdef MPI
/** Since datasets are actually # frames * however many atoms found in mask
  * at each frame, sync here.
  */
int Action_Mask::SyncAction(Parallel::Comm const& commIn) {
  if (fnum_ == 0) return 0;
  // Get total number of mask entries.
  std::vector<int> rank_frames( commIn.Size() );
  commIn.GatherMaster( &idx_, 1, MPI_INT, &(rank_frames[0]) );
  mprintf("DEBUG: Master= %i frames\n", idx_);
  for (int rank = 1; rank < commIn.Size(); rank++) {
    mprintf("DEBUG: Rank%i= %i frames\n", rank, rank_frames[ rank ]);
    idx_ += rank_frames[ rank ];
  }
  mprintf("DEBUG: Total= %i frames.\n", idx_);
  fnum_->Sync( idx_, rank_frames, commIn );
  fnum_->SetNeedsSync( false );
  anum_->Sync( idx_, rank_frames, commIn );
  anum_->SetNeedsSync( false );
  aname_->Sync( idx_, rank_frames, commIn );
  aname_->SetNeedsSync( false );
  rnum_->Sync( idx_, rank_frames, commIn );
  rnum_->SetNeedsSync( false );
  rname_->Sync( idx_, rank_frames, commIn );
  rname_->SetNeedsSync( false );
  mnum_->Sync( idx_, rank_frames, commIn );
  mnum_->SetNeedsSync( false );

  return 0;
}
#endif
