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
  axis_(DZ),
  //area_coord_{DZ, DY},		// this is C++11!
  property_(NUMBER),
  delta_(0.0)
{
  area_coord_[0] = DX; area_coord_[1] = DY;
}

void Action_Density::Help()
{
  mprintf("\tout <filename> [delta <resolution>] [x|y|z]\n"
	  "\t[number|mass|charge|electron] <mask1> ... <maskN>\n"
          "  Calculate density along a coordinate.\n");
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

  if (output_.OpenEnsembleWrite(outfileName, DSL->EnsembleNum()) ) {
    mprinterr("Error: Density: Could not open output file %s\n",
	      outfileName.c_str());
    return Action::ERR;
  }

  if (actionArgs.hasKey("x") ) {
    axis_ = DX;
    area_coord_[0] = DY;
    area_coord_[1] = DZ;
  } else if (actionArgs.hasKey("y") ) {
    axis_ = DY;
    area_coord_[0] = DX;
    area_coord_[1] = DZ;
  } else if (actionArgs.hasKey("z") ) {
    axis_ = DZ;
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
    masks_.push_back( AtomMask(maskstr) );
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
      slice = (long) (coord[axis_] / delta_);

      // FIXME: split 0 bin in one + and one -
      if (slice != 0)
	tmp[slice] += properties_[i][j];
      else
	if (coord[axis_] < 0.0)
	  tmp[slice-1] += properties_[i][j];
	else
	  tmp[slice+1] += properties_[i][j];

      j++;
    }

    if (tmp.size() > 0)
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

  bool first, scale_area;
  long minidx = LONG_MAX, maxidx = LONG_MIN;
  double density, sd, area;

  std::map<int,double>::iterator e, b, itv;
  statmap curr;



  area = area_.mean();
  sd = sqrt(area_.variance());
  scale_area = (property_ == ELECTRON && area > SMALL);

  mprintf("The average box area in %c/%c is %.2f Angstrom (sd = %.2f).\n",
	  area_coord_[0] + 88, area_coord_[1] + 88, area, sd);

  if (scale_area)
    mprintf("The electron density will be scaled by this area.\n");

  // the search for minimum and maximum indices relies on ordered map
  for (unsigned long i = 0; i < histograms_.size(); i++) {
    b = histograms_[i].mean_begin(); 
    e = histograms_[i].mean_end();

    if (b->first < minidx)
      minidx = b->first;

    if (e != b) {
      e--;
      if (e->first > maxidx)
	maxidx = e->first;
    }
  }

  output_.Printf("#density");

  for (std::vector<AtomMask>::const_iterator mask = masks_.begin();
       mask != masks_.end();
       mask++) {
    output_.Printf(" %s sd(%s)", mask->MaskString(), mask->MaskString() );
  }

  output_.Printf("\n");

  // make sure we have zero values at beginning and end as this
  // "correctly" integrates the histogram
  minidx--;
  maxidx++;

  for (long i = minidx; i <= maxidx; i++) {
    first = true;
    if (i == 0) continue;	// FIXME: 0 is doubly counted

    for (unsigned long j = 0; j < histograms_.size(); j++) {
      curr = histograms_[j];

      if (first) {
        output_.Printf("%10.4f", (i < 0 ? -delta_ : 0.0) +
                       ((double) i + 0.5) * delta_);
        first = false;
      }

      density = curr.mean(i) / (delta_ * // FIXME: 0 is doubly counted
				(i == -1 or i == 1 ? 2.0 : 1.0 ));
      sd = sqrt(curr.variance(i) );

      if (scale_area) {
	density /= area;
	sd /= area;
      }

      output_.Printf(" %10.3f %10.5f", density, sd);
    }

    output_.Printf("\n");
  }
}

