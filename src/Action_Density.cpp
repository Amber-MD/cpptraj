#include <cmath>
#include <climits>

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
  //area_coord_{DZ, DY},		// this is C++11!
  property_(NUMBER),
  delta_(0.0)
{
  area_coord_[0] = DX; area_coord_[1] = DY;
}

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

  if (actionArgs.hasKey("x") ) {
    direction_ = DX;
    area_coord_[0] = DY;
    area_coord_[1] = DZ;
  } else if (actionArgs.hasKey("y") ) {
    direction_ = DY;
    area_coord_[0] = DX;
    area_coord_[1] = DZ;
  } else if (actionArgs.hasKey("z") ) {
    direction_ = DZ;
    area_coord_[0] = DX;
    area_coord_[1] = DY;
  }

  property_ = NUMBER;
  if (actionArgs.hasKey("number") )   property_ = NUMBER;
  if (actionArgs.hasKey("mass") )     property_ = MASS;
  if (actionArgs.hasKey("charge") )   property_ = CHARGE;
  if (actionArgs.hasKey("electron") ) property_ = ELECTRON;

  delta_ = actionArgs.getKeyDouble("delta", 0.01);

  // for compatibility with ptraj, ignored because we rely on the atom code to
  // do the right thing, see Atom.{h,cpp}
 actionArgs.GetStringKey("efile");

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
      case NUMBER:
	property.push_back(1.0);
	break;

      case MASS:
	property.push_back(atom.Mass() );
	break;

      case CHARGE:
	property.push_back(atom.Charge() );
	break;

      case ELECTRON:
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
  long slice;
  unsigned long i, j;
  Vec3 coord;
  Box box;


  i = 0;

  for (std::vector<AtomMask>::const_iterator mask = masks_.begin();
       mask != masks_.end();
       mask++) {

    j = 0;

    std::map<int,double> tmp;

    for (AtomMask::const_iterator idx = mask->begin();
	 idx != mask->end();
	 idx++) {
      coord = currentFrame->XYZ(*idx);
      slice = (long) (coord[direction_] / delta_);

      // FIXME: split 0 bin in one + and one -
      tmp[slice] += properties_[i][j];

      j++;
    }

    histograms_[i].accumulate(tmp);
    i++;
  }

  box = currentFrame->BoxCrd();
  area_.accumulate(box[area_coord_[0]] * box[area_coord_[1]]);

  return Action::OK;
}


// Action_Density::print()
void Action_Density::Print()
{
  const unsigned int SMALL = 1.0;

  bool first;
  long minidx = LONG_MAX, maxidx = LONG_MIN;
  double density, sd, area;

  std::map<int,double>::iterator it, itv;
  statmap curr;



  area = area_.mean();

  mprintf("The average box area in %c/%c is %.2f Angstrom (sd = %.2f).\n",
	  area_coord_[0] + 88, area_coord_[1] + 88, area,
	  sqrt(area_.variance()) );

  if (property_ == ELECTRON && area > SMALL)
    mprintf("The electron density will be scaled by this area.\n");

  // the search for minimum and maximum indices relies on ordered map
  for (unsigned long i = 0; i < histograms_.size(); i++) {
    it = histograms_[i].mean_begin(); 
    if (it->first < minidx)
      minidx = it->first;

    it = histograms_[i].mean_end();
    it--;
    if (it->first > maxidx)
      maxidx = it->first;
  }

  output_.Printf("#dist");

  for (std::vector<AtomMask>::const_iterator mask = masks_.begin();
       mask != masks_.end();
       mask++) {
    output_.Printf(" %s sd(%s)", mask->MaskString(), mask->MaskString() );
  }

  output_.Printf("\n");

  for (long i = minidx; i <= maxidx; i++) {
    first = true;

    for (unsigned long j = 0; j < histograms_.size(); j++) {
      curr = histograms_[j];

      if (first) {
        output_.Printf("%10.4f", (i < 0 ? -delta_ : 0.0) +
                       ((double) i + (i == 0 ? 0.0 : 0.5)) * delta_);
        first = false;
      }

      density = curr.mean(i) / (delta_ * // 0 is double counted
				(i == 0 ? 2.0 : 1.0 ));
      sd = sqrt(curr.variance(i) );

      if (property_ == ELECTRON && area > SMALL) {
	density /= area;
	sd /= area;
      }

      output_.Printf(" %10.3f %10.5f", density, sd);
    }

    output_.Printf("\n");
  }
}

