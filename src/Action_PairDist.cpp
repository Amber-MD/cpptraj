#include <cmath>
#include "Action_PairDist.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"
#include "DistRoutines.h"


/** Calculate pair distribution function P(r) between two masks.
  * \author Hannes H. Loeffler.
  */

// CONSTRUCTOR
Action_PairDist::Action_PairDist() :
  Pr_(0),
  delta_(0.0),
  maxbin_(-1)
{}

void Action_PairDist::Help()
{
  mprintf("\tout <filename> mask <mask> [mask2 <mask>] [<resolution> delta]\n");
}

// Action_PairDist::init()
Action::RetType Action_PairDist::Init(ArgList& actionArgs,
				      TopologyList* PFL, FrameList* FL,
				      DataSetList* DSL, DataFileList* DFL,
				      int debugIn)
{
  InitImaging(true);

  std::string outfileName = actionArgs.GetStringKey("out");

  if (outfileName.empty()) {
    outfileName = "pairdist.dat";
  }

  std::string mask1 = actionArgs.GetStringKey("mask");

  if (mask1.empty()) {
    mprinterr("Error: pairdist: No mask1 specified.\n");
    return Action::ERR;
  }

  mask1_.SetMaskString(mask1);

  std::string mask2 = actionArgs.GetStringKey("mask2");

  if (mask2.empty()) {
    mask2_.SetMaskString(mask1);
  } else {
    mask2_.SetMaskString(mask2);
  }

  delta_ = actionArgs.getKeyDouble("delta", 0.01);

  Pr_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), "P(r)");
  Pr_->SetPrecision(12, 0);
  DataFile* outfile = DFL->AddSetToFile(outfileName, Pr_);

  if (outfile == 0) {
    mprinterr("Error: PairDist: Could not setup output file %s\n",
	      outfileName.c_str());
    return Action::ERR;
  }

  outfile->ProcessArgs("xmin " + doubleToString(delta_ / 2.0) + 
		       " xstep " + doubleToString(delta_) );

  return Action::OK;
}


// Action_PairDist::Setup()
Action::RetType Action_PairDist::Setup(Topology* currentParm,
				       Topology** parmAddress)
{
  if (currentParm->SetupIntegerMask(mask1_) ) return Action::ERR;

  mprintf("\t");
  mask1_.BriefMaskInfo();

  if (mask1_.None()) {
    mprintf("    Error: PairDist::setup: Mask has no atoms.\n");
    return Action::ERR;
  }

  if (currentParm->SetupIntegerMask(mask2_) ) return Action::ERR;

  mask2_.BriefMaskInfo();
  mprintf("\n");

  if (mask2_.None()) {
    mprintf("    Error: PairDist::setup: Mask2 has no atoms.\n");
    return Action::ERR;
  }

  SetupImaging(currentParm->BoxType() );

  return Action::OK;  
}


// Action_PairDist::action()
Action::RetType Action_PairDist::DoAction(int frameNum,
					  Frame* currentFrame,
					  Frame** frameAddress)
{
  int bin;
  double Dist = 0.0;
  Matrix_3x3 ucell, recip;
  Vec3 a1, a2;


  for (AtomMask::const_iterator atom1 = mask1_.begin();
       atom1 != mask1_.end(); atom1++) {
    for (AtomMask::const_iterator atom2 = mask2_.begin();
	 atom2 != mask2_.end(); atom2++) {

      if (*atom1 == *atom2) continue;

      a1 = currentFrame->XYZ(*atom1);
      a2 = currentFrame->XYZ(*atom2);

      switch (ImageType() ) {
      case NONORTHO:
	currentFrame->BoxCrd().ToRecip(ucell, recip);
	Dist = DIST2_ImageNonOrtho(a1, a2, ucell, recip);
	break;
      case ORTHO:
	Dist = DIST2_ImageOrtho(a1, a2, currentFrame->BoxCrd());
	break;
      case NOIMAGE:
	Dist = DIST2_NoImage(a1, a2);
	break;
      }

      bin = (int) (sqrt(Dist) / delta_);

      if (bin > maxbin_) {
	maxbin_ = bin;
	histogram_.resize(maxbin_ + 1);
      }

      histogram_[bin] += 1.0;
    }
  }

  return Action::OK;
}


// Action_PairDist::print()
void Action_PairDist::Print()
{
  int i;
  std::vector<double>::iterator it;


  for (i = 0, it = histogram_.begin();
       it != histogram_.end(); ++it) {
    Pr_->Add(i++, &(*it) );
  }
}
