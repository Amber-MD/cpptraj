#include <cmath>

#include "Action_Density.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"
#include "DistRoutines.h"
#include "DataSet_Mesh.h"


/** Calculate density along a coordinate.
  * \author Hannes H. Loeffler.
  */

const std::string Action_Density::emptystring = "";

// CONSTRUCTOR
Action_Density::Action_Density() :
  axis_(DZ),
  //area_coord_{DZ, DY},    // this is C++11!
  property_(NUMBER),
  delta_(0.0)
{
  area_coord_[0] = DX; area_coord_[1] = DY;
}

void Action_Density::Help() const {
  mprintf("\t[out <filename>] [name <set name>] [delta <resolution>] [x|y|z]\n"
          "\t[number|mass|charge|electron] <mask1> ... <maskN>\n"
          "  Calculate density along a coordinate.\n");
}

const char* Action_Density::PropertyStr_[] = {
  "number", "mass", "charge", "electron"
};

const char Action_Density::AxisStr_[] = { 'X', 'Y', 'Z' };

// Action_Density::Init()
Action::RetType Action_Density::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  if (init.TrajComm().Size() > 1) {
    mprinterr("Error: 'density' action does not work with > 1 thread (%i threads currently).\n",
              init.TrajComm().Size());
    return Action::ERR;
  }
# endif
  InitImaging(true);

  DataFile* outfile = init.DFL().AddDataFile(actionArgs.GetStringKey("out"), actionArgs);

  std::string dsname = actionArgs.GetStringKey("name");
  if (dsname.empty())
    dsname = init.DSL().GenerateDefaultName("DENSITY");

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
  if (actionArgs.hasKey("efile"))
    mprintf("Warning: The 'efile' keyword is deprecated.\n");

  // read the rest of the command line as a series of masks
  std::string maskstr;

  unsigned int idx = 1;
  while ( (maskstr = actionArgs.GetMaskNext() ) != emptystring) {
    masks_.push_back( AtomMask(maskstr) );
    MetaData MD(dsname, "avg", idx);
    MD.SetTimeSeries( MetaData::NOT_TS );
    // Hold average density
    DataSet* ads = init.DSL().AddSet( DataSet::XYMESH, MD );
    if (ads == 0) return Action::ERR;
    ads->SetLegend( NoWhitespace(masks_.back().MaskExpression()) );
    AvSets_.push_back( ads );
    if (outfile != 0) outfile->AddDataSet( ads );
    // Hold SD density
    MD.SetAspect("sd");
    DataSet* sds = init.DSL().AddSet( DataSet::XYMESH, MD );
    if (sds == 0) return Action::ERR;
    sds->SetLegend( NoWhitespace("sd(" + masks_.back().MaskExpression() + ")") );
    SdSets_.push_back( sds );
    if (outfile != 0) outfile->AddDataSet( sds );
#   ifdef MPI
    ads->SetNeedsSync( false ); // Populated in Print()
    sds->SetNeedsSync( false );
#   endif
    idx++;
  }
  if (masks_.empty()) {
    mprinterr("Error: No masks specified.\n");
    return Action::ERR;
  }

  minus_histograms_.resize(masks_.size() );
  plus_histograms_.resize(masks_.size() );

  mprintf("    DENSITY: Determining %s density for %zu masks.\n", PropertyStr_[property_],
          masks_.size());
  mprintf("\troutine version: %s\n", ROUTINE_VERSION_STRING);
  mprintf("\tDelta is %f\n", delta_);
  mprintf("\tAxis is %c\n", AxisStr_[axis_]);
  mprintf("\tData set name is '%s'\n", dsname.c_str());
  if (outfile != 0)
    mprintf("\tOutput to '%s'\n", outfile->DataFilename().full());

  return Action::OK;
}

// Action_Density::Setup()
Action::RetType Action_Density::Setup(ActionSetup& setup) {
  properties_.resize(0);

  for (std::vector<AtomMask>::iterator mask = masks_.begin();
       mask != masks_.end();
       mask++) {
    if (setup.Top().SetupIntegerMask(*mask) ) return Action::ERR;

    std::vector<double> property;

    for (AtomMask::const_iterator idx = mask->begin();
      idx != mask->end(); idx++)
    {
      const Atom& atom = setup.Top()[*idx];

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

  SetupImaging(setup.CoordInfo().TrajBox().Type() );

  return Action::OK;  
}

// Action_Density::DoAction()
Action::RetType Action_Density::DoAction(int frameNum, ActionFrame& frm) {
  long slice;
  unsigned long i, j;
  Vec3 coord;
  Box box;

  i = 0;

  for (std::vector<AtomMask>::const_iterator mask = masks_.begin();
       mask != masks_.end();
       mask++)
  {
    j = 0;

    std::map<long,double> minus_histo, plus_histo;

    for (AtomMask::const_iterator idx = mask->begin();
         idx != mask->end();
         idx++)
    {
      coord = frm.Frm().XYZ(*idx);
      slice = (unsigned long) (coord[axis_] / delta_);

      if (coord[axis_] < 0) {
        minus_histo[slice] += properties_[i][j];
      } else {
        plus_histo[slice] += properties_[i][j];
      }

      j++;
    }

    if (minus_histo.size() > 0)
      minus_histograms_[i].accumulate(minus_histo);

    if (plus_histo.size() > 0)
      plus_histograms_[i].accumulate(plus_histo);

    i++;
  }

  box = frm.Frm().BoxCrd();
  area_.accumulate(box[area_coord_[0]] * box[area_coord_[1]]);

  return Action::OK;
}


// Action_Density::Print()
void Action_Density::Print()
{
  const unsigned int SMALL = 1.0;

  //bool first_round, scale_area;
  long minus_minidx = 0, minus_maxidx = 0, plus_minidx = 0, plus_maxidx = 0;
  double density, sd, area;

  std::map<long,double>::iterator first_idx, last_idx;
  statmap curr;



  area = area_.mean();
  sd = sqrt(area_.variance());
  bool scale_area = (property_ == ELECTRON && area > SMALL);

  mprintf("    DENSITY: The average box area in %c/%c is %.2f Angstrom (sd = %.2f).\n",
          area_coord_[0] + 88, area_coord_[1] + 88, area, sd);

  if (scale_area)
    mprintf("The electron density will be scaled by this area.\n");

  // the search for minimum and maximum indices relies on ordered map
  for (unsigned long i = 0; i < minus_histograms_.size(); i++) {
    first_idx = minus_histograms_[i].mean_begin(); 
    last_idx = minus_histograms_[i].mean_end();

    if (first_idx->first < minus_minidx)
      minus_minidx = first_idx->first;

    if (last_idx != first_idx) {
      last_idx--;
      if (last_idx->first > minus_maxidx)
        minus_maxidx = last_idx->first;
    }
  }

  for (unsigned long i = 0; i < plus_histograms_.size(); i++) {
    first_idx = plus_histograms_[i].mean_begin(); 
    last_idx = plus_histograms_[i].mean_end();

    if (first_idx->first < plus_minidx)
      plus_minidx = first_idx->first;

    if (last_idx != first_idx) {
      last_idx--;
      if (last_idx->first > plus_maxidx)
        plus_maxidx = last_idx->first;
    }
  }

  //output_->Printf("# routine version: %s\n#density", ROUTINE_VERSION_STRING);

  //for (std::vector<AtomMask>::const_iterator mask = masks_.begin();
  //     mask != masks_.end();
  //     mask++) {
  //  output_->Printf(" %s sd(%s)", mask->MaskString(), mask->MaskString() );
  //}

  //output_->Printf("\n");

  // make sure we have zero values at beginning and end as this
  // "correctly" integrates the histogram
  minus_minidx--;
  plus_maxidx++;

  // Set up data set dimensions
  double Xmin = -delta_ + ((double) minus_minidx + 0.5) * delta_;
  Dimension Xdim(Xmin, delta_, "density");
  for (unsigned int j = 0; j != AvSets_.size(); j++) {
    AvSets_[j]->SetDim(Dimension::X, Xdim);
    SdSets_[j]->SetDim(Dimension::X, Xdim);
  }

  for (long i = minus_minidx; i <= minus_maxidx; i++) {
    //first_round = true;
    double Xcoord = -delta_ + ((double) i + 0.5) * delta_;

    for (unsigned long j = 0; j < minus_histograms_.size(); j++) {
      curr = minus_histograms_[j];

      //if (first_round) {
      //  output_->Printf("%10.4f", -delta_ + ((double) i + 0.5) * delta_);
      //  first_round = false;
      //}
      
      density = curr.mean(i) / delta_;
      sd = sqrt(curr.variance(i) );

      if (scale_area) {
        density /= area;
        sd /= area;
      }

      //output_->Printf(" %10.3f %10.5f", density, sd);
      ((DataSet_Mesh*)AvSets_[j])->AddXY( Xcoord, density );
      ((DataSet_Mesh*)SdSets_[j])->AddXY( Xcoord, sd      );
    }

    //output_->Printf("\n");
  }

  for (long i = plus_minidx; i <= plus_maxidx; i++) {
    //first_round = true;
    double Xcoord = ((double) i + 0.5) * delta_;

    for (unsigned long j = 0; j < plus_histograms_.size(); j++) {
      curr = plus_histograms_[j];

      //if (first_round) {
      //  output_->Printf("%10.4f", ((double) i + 0.5) * delta_);
      //  first_round = false;
      //}

      density = curr.mean(i) / delta_;
      sd = sqrt(curr.variance(i) );

      if (scale_area) {
        density /= area;
        sd /= area;
      }

      //output_->Printf(" %10.3f %10.5f", density, sd);
      ((DataSet_Mesh*)AvSets_[j])->AddXY( Xcoord, density );
      ((DataSet_Mesh*)SdSets_[j])->AddXY( Xcoord, sd      );
    }

    //output_->Printf("\n");
  }
}
