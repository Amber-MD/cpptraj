// Action_MinDist
#include <cmath>
#include <cfloat> // DBL_MAX
#include "Action_MinDist.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_MinDist::Action_MinDist() :
  dist_(0)
{ } 

void Action_MinDist::Help() {
  mprintf("\t[<name>] <mask1> <mask2> [out <filename>] [noimage]\n");
}

// Action_MinDist::Init()
Action::RetType Action_MinDist::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get Keywords
  image_.InitImaging( !(actionArgs.hasKey("noimage")) );
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );

  // Get Masks
  std::string mask1 = actionArgs.GetMaskNext();
  std::string mask2 = actionArgs.GetMaskNext();
  if (mask1.empty() || mask2.empty()) {
    mprinterr("Error: mindist: Requires 2 masks\n");
    return Action::ERR;
  }
  Mask1_.SetMaskString(mask1);
  Mask2_.SetMaskString(mask2);

  // Dataset to store minimum distances
  dist_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), "MinDist");
  if (dist_==0) return Action::ERR;
  // Add dataset to data file
  if (outfile != 0) outfile->AddSet( dist_ );

  mprintf("    MINDIST: looking for minimum distance of atoms in %s to atoms in %s\n\t",
           Mask1_.MaskString(), Mask2_.MaskString());
  if (!image_.UseImage()) 
    mprintf("Non-imaged");
  else
    mprintf("Imaged");
  mprintf(".\n");

  return Action::OK;
}

// Action_MinDist::Setup()
/** Determine what atoms each mask pertains to for the current parm file.
  */
Action::RetType Action_MinDist::Setup(Topology* currentParm, Topology** parmAddress) {
  if (currentParm->SetupIntegerMask( Mask1_ )) return Action::ERR;
  if (currentParm->SetupIntegerMask( Mask2_ )) return Action::ERR;
  Mask1_.MaskInfo();
  Mask2_.MaskInfo();
  if (Mask1_.None() || Mask2_.None()) {
    mprintf("Warning: mindist: One or both masks have no atoms.\n");
    return Action::ERR;
  }
  // Check if the masks overlap at all
  if ( Mask1_.NumAtomsInCommon( Mask2_ ) > 0 ) {
    mprintf("Error: mindist: Masks overlap, minimum distance will always be 0.0\n");
    return Action::ERR;
  }
  // Set up imaging info for this parm
  image_.SetupImaging( currentParm->BoxType() );
  if (image_.ImagingEnabled())
    mprintf("\tImaged.\n");
  else
    mprintf("\tImaging off.\n");

  return Action::OK;  
}

// Action_MinDist::DoAction()
Action::RetType Action_MinDist::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  double Dist;
  AtomMask::const_iterator atom1, atom2;
  Matrix_3x3 ucell, recip;
  double mindist2 = DBL_MAX;

  switch ( image_.ImageType() ) {
    case NONORTHO:
      currentFrame->BoxCrd().ToRecip(ucell, recip);
      for (atom1 = Mask1_.begin(); atom1 != Mask1_.end(); ++atom1) {
        for (atom2 = Mask2_.begin(); atom2 != Mask2_.end(); ++atom2) {
          Dist = DIST2_ImageNonOrtho(currentFrame->XYZ(*atom1), currentFrame->XYZ(*atom2), 
                                     ucell, recip);
          if (Dist < mindist2) mindist2 = Dist;
        }
      }
      break;
    case ORTHO:
      for (atom1 = Mask1_.begin(); atom1 != Mask1_.end(); ++atom1) {
        for (atom2 = Mask2_.begin(); atom2 != Mask2_.end(); ++atom2) {
          Dist = DIST2_ImageOrtho(currentFrame->XYZ(*atom1), currentFrame->XYZ(*atom2),
                                  currentFrame->BoxCrd());
          if (Dist < mindist2) mindist2 = Dist;
        }
      }
      break;
    case NOIMAGE:
      for (atom1 = Mask1_.begin(); atom1 != Mask1_.end(); ++atom1) {
        for (atom2 = Mask2_.begin(); atom2 != Mask2_.end(); ++atom2) {
          Dist = DIST2_NoImage(currentFrame->XYZ(*atom1), currentFrame->XYZ(*atom2));
          if (Dist < mindist2) mindist2 = Dist;
        }
      }
      break;
  }
  Dist = sqrt(mindist2);

  dist_->Add(frameNum, &Dist);

  //fprintf(outfile,"%10i %10.4lf\n",frameNum,D);
  
  return Action::OK;
}
