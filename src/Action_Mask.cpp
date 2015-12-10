// Action_Mask
#include "Action_Mask.h"
#include "CpptrajStdio.h"
#include "Trajout_Single.h"

// CONSTRUCTOR
Action_Mask::Action_Mask() :
  ensembleNum_(-1),
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

  mprintf("    ACTIONMASK: Information on atoms in mask %s will be printed",
          Mask1_.MaskString());
  if (outfile_ != 0)
    mprintf(" to file %s",outfile_->Filename().full());
  mprintf(".\n");
  if (!maskpdb_.empty()) 
    mprintf("\t%ss of atoms in mask will be written to %s.X\n",
            TrajectoryFile::FormatString(trajFmt_), maskpdb_.c_str());
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
      if (outfile_ != 0)
        outfile_->Printf("%8i %8i %4s %8i %4s %8i\n", frameNum+1,
                        atom+1, (*CurrentParm_)[atom].c_str(), res+1,
                        CurrentParm_->Res(res).c_str(), (*CurrentParm_)[atom].MolNum()+1);
      /*mprintf(" Type=%4s",CurrentParm_->types[atom]);
      mprintf(" Charge=%lf",CurrentParm_->charge[atom]);
      mprintf(" Mass=%lf",CurrentParm_->mass[atom]);
      outfile.Printf("\n");*/
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
