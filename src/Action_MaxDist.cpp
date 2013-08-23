// Action_MaxDist
#include <cmath>
#include "Action_MaxDist.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_MaxDist::Action_MaxDist() :
  dist_(0)
{ } 

void Action_MaxDist::Help() {
  mprintf("\t[<name>] <mask1> [<mask2>] [out <filename>] [noimage]\n");
}

// Action_MaxDist::Init()
Action::RetType Action_MaxDist::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get Keywords
  image_.InitImaging( !(actionArgs.hasKey("noimage")) );
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );

  // Get Masks
  std::string mask1 = actionArgs.GetMaskNext();
  if (mask1.empty()) {
    mprinterr("Error: maxdist: Requires at least 1 mask\n");
    return Action::ERR;
  }
  Mask1_.SetMaskString(mask1);
  std::string mask2 = actionArgs.GetMaskNext();
  if (!mask2.empty()) Mask2_.SetMaskString(mask2);

  // Dataset to store max distances
  dist_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), "MaxDist");
  if (dist_==0) return Action::ERR;
  // Add dataset to data file
  if (outfile != 0) outfile->AddSet( dist_ );

  mprintf("    MAXDIST: looking for max distance");
  if (Mask2_.MaskStringSet())
    mprintf(" of atoms in %s to atoms in %s\n", Mask1_.MaskString(), Mask2_.MaskString());
  else
    mprintf(" between atoms in %s\n", Mask1_.MaskString());
  if (!image_.UseImage()) 
    mprintf("\tNon-imaged.\n");
  else
    mprintf("\tImaged.\n");

  return Action::OK;
}

// Action_MaxDist::Setup()
/** Determine what atoms each mask pertains to for the current parm file.
  */
Action::RetType Action_MaxDist::Setup(Topology* currentParm, Topology** parmAddress) {
  if (currentParm->SetupIntegerMask( Mask1_ )) return Action::ERR;
  Mask1_.MaskInfo();
  if (Mask1_.None()) {
    mprintf("Warning: maxdist: First mask has no atoms.\n");
    return Action::ERR;
  }
  if (Mask2_.MaskStringSet()) {
    if (currentParm->SetupIntegerMask( Mask2_ )) return Action::ERR;
    Mask2_.MaskInfo();
    if (Mask2_.None()) {
      mprintf("Warning: maxdist: Second mask has no atoms.\n");
      return Action::ERR;
    }
    // Warn if the masks overlap at all
    if ( Mask1_.NumAtomsInCommon( Mask2_ ) > 0 )
      mprintf("Warning: maxdist: Masks overlap.\n");
  }
  // Set up imaging info for this parm
  image_.SetupImaging( currentParm->BoxType() );
  if (image_.ImagingEnabled())
    mprintf("\tImaged.\n");
  else
    mprintf("\tImaging off.\n");

  return Action::OK;  
}

// Action_MaxDist::DoAction()
Action::RetType Action_MaxDist::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  Matrix_3x3 ucell, recip;
  double maxdist2 = 0.0;

  if (image_.ImageType() == NONORTHO) currentFrame->BoxCrd().ToRecip(ucell, recip);
  if (Mask2_.MaskStringSet()) {
    
    for (AtomMask::const_iterator atom1 = Mask1_.begin(); atom1 != Mask1_.end(); ++atom1) {
      for (AtomMask::const_iterator atom2 = Mask2_.begin(); atom2 != Mask2_.end(); ++atom2) {
        double Dist = DIST2(currentFrame->XYZ(*atom1), currentFrame->XYZ(*atom2),
                            image_.ImageType(), currentFrame->BoxCrd(), ucell, recip);
        if (Dist > maxdist2) maxdist2 = Dist;
      }
    }
  } else {
    AtomMask::const_iterator atom1_end = Mask1_.end() - 1;
    for (AtomMask::const_iterator atom1 = Mask1_.begin(); atom1 != atom1_end; ++atom1) {
      for (AtomMask::const_iterator atom2 = atom1 + 1; atom2 != Mask1_.end(); ++atom2) {
        double Dist = DIST2(currentFrame->XYZ(*atom1), currentFrame->XYZ(*atom2),
                            image_.ImageType(), currentFrame->BoxCrd(), ucell, recip);
        if (Dist > maxdist2) maxdist2 = Dist;
      }
    }
  }
  maxdist2 = sqrt(maxdist2);

  dist_->Add(frameNum, &maxdist2);

  //fprintf(outfile,"%10i %10.4lf\n",frameNum,D);
  
  return Action::OK;
}
