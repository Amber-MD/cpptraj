#include <cmath>

#include "Action_Density.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"
#include "DistRoutines.h"


/** Calculate density along a coordinate.
  * \author Hannes H. Loeffler.
  */

const std::string Action_Density::emptystring = "";

// CONSTRUCTOR
Action_Density::Action_Density() :
  direction_(DZ),
  property_(NUMBER),
  frameNum_(0),
  delta_(0.0)
{}

void Action_Density::Help()
{
  mprintf("\tout <filename> [delta <resolution>] [x|y|z]\n"
	  "[number|mass|charge|electron] <mask1> ... <maskN>\n");
}

// Action_Density::init()
Action::RetType Action_Density::Init(ArgList& actionArgs,
				     TopologyList* PFL, FrameList* FL,
				     DataSetList* DSL, DataFileList* DFL,
				     int debugIn)
{
  InitImaging(true);

  std::string outfileName = actionArgs.GetStringKey("out");

  if (outfileName.empty()) {
    outfileName = "density.dat";
  }

  if (output_.OpenWrite(outfileName) ) {
    mprinterr("Error: Density: Could not open output file %s\n",
	      outfileName.c_str());
    return Action::ERR;
  }

  direction_ = DZ;
  if (actionArgs.hasKey("x") ) direction_ = DX;
  if (actionArgs.hasKey("y") ) direction_ = DY;
  if (actionArgs.hasKey("z") ) direction_ = DZ;

  property_ = NUMBER;
  if (actionArgs.hasKey("number") )   property_ = NUMBER;
  if (actionArgs.hasKey("mass") )     property_ = MASS;
  if (actionArgs.hasKey("charge") )   property_ = CHARGE;
  if (actionArgs.hasKey("electron") ) property_ = ELECTRON;

  delta_ = actionArgs.getKeyDouble("delta", 0.01);

  // read the rest of the command line as a series of masks
  std::string maskstr;

  while ( (maskstr = actionArgs.GetMaskNext() ) != emptystring) {
    masks_.push_back(AtomMask(maskstr) );
  }

  histograms_.resize(masks_.size() );

  return Action::OK;
}


// Action_Density::Setup()
Action::RetType Action_Density::Setup(Topology* currentParm,
				      Topology** parmAddress)
{
  properties_.resize(0);

  for (std::vector<AtomMask>::iterator mask = masks_.begin();
       mask != masks_.end();
       mask++) {
    if (currentParm->SetupIntegerMask(*mask) ) return Action::ERR;

    std::vector<double> property;

    for (AtomMask::const_iterator idx = mask->begin();
	 idx != mask->end(); idx++) {
      const Atom& atom = (*currentParm)[*idx];

      switch (property_) {
      case NUMBER:		// FIXME: probably redundant, check in DoAction
	property.push_back(1.0);
	break;

      case MASS:
	property.push_back(atom.Mass() );
	break;

      case CHARGE:
	property.push_back(atom.Charge() );
	break;

      case ELECTRON:
	// rely on the presence of AtomicNumber, ptraj reads efile
	property.push_back(atom.AtomicNumber() - atom.Charge() );
	break;
      }
    }

    properties_.push_back(property);

    mprintf("\t");
    mask->BriefMaskInfo();
    mprintf("\n");
  }

  SetupImaging(currentParm->BoxType() );

  return Action::OK;  
}


// Action_Density::action()
Action::RetType Action_Density::DoAction(int frameNum,
					  Frame* currentFrame,
					  Frame** frameAddress)
{
  int slice, i, j;
  Vec3 coord;


  i = 0;

  for (std::vector<AtomMask>::iterator mask = masks_.begin();
       mask != masks_.end();
       mask++) {
    std::map<int,double> h = histograms_[i];

    j = 0;

    for (AtomMask::const_iterator idx = mask->begin();
	 idx != mask->end(); idx++) {
      coord = currentFrame->XYZ(*idx);
      slice = (int) (coord[direction_] / delta_);

      h[slice] += properties_[i][j];

      j++;
    }

    histograms_[i] = h;

    i++;
  }

  frameNum_++;			// FIXME: check if otherwise available

  return Action::OK;
}


// Action_Density::print()
void Action_Density::Print()
{
  std::map<int,double>::iterator it;


  for (unsigned int i = 0; i < histograms_.size(); i++) {
    output_.Printf("#%s\n", masks_[i].MaskString() );

    for (it = histograms_[i].begin();
	 it != histograms_[i].end(); ++it) {
      output_.Printf("%10.4f %10.3f\n",
		     // FIXME: check if scaling correct, currently same as ptraj
		     (it->first < 0 ? -delta_ : 0.0) +
		     ((double) it->first + 0.5) * delta_,
		     it->second / (delta_ * frameNum_ * // 0 is double counted
				   (it->first == 0 ? 2.0 : 1.0 )) );
    }

    output_.Printf("\n");
  }
}
