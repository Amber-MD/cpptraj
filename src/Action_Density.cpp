#include <algorithm> // std::min, std::max
#include <cmath>
#include "Action_Density.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // NoWhitespace()
#include "Constants.h"
#include "DataSet_1D.h"

const std::string Action_Density::emptystring = "";

// CONSTRUCTOR
Action_Density::Action_Density() :
  axis_(DZ),
  property_(NUMBER),
  binType_(CENTER),
  delta_(0.0),
  density_(0),
  restrictType_(NONE),
  cutVal_(-1)
{
  area_coord_[0] = DX; area_coord_[1] = DY;
}

void Action_Density::Help() const {
  mprintf("\t[out <filename>] [name <set name>]\n"
          "\t[ <mask1> ... <maskN> [delta <resolution>] [{x|y|z}]\n"
          "\t  [{number|mass|charge|electron}] [{bincenter|binedge}]\n"
          "\t  [restrict {cylinder|square} cutoff <cut>]  ]\n"
          "  If one or more masks are specified, calculate specified density of selected\n"
          "  atoms along a coordinate. Otherwise calculate the total system density.\n");
}

const char* Action_Density::PropertyStr_[] = {
  "number", "mass", "charge", "electron"
};

const char* Action_Density::AxisStr_[] = { "X", "Y", "Z" };

// Action_Density::Init()
Action::RetType Action_Density::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  trajComm_ = init.TrajComm();
# endif
  DataFile* outfile = init.DFL().AddDataFile(actionArgs.GetStringKey("out"), actionArgs);

  std::string dsname = actionArgs.GetStringKey("name");

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

  restrictType_ = NONE;
  std::string restrictStr = actionArgs.GetStringKey("restrict");
  if (!restrictStr.empty()) {
    if (restrictStr == "cylinder")
      restrictType_ = CYLINDER;
    else if (restrictStr == "square")
      restrictType_ = SQUARE;
    else {
      mprinterr("Error: Unrecognized option for 'restrict': %s\n", restrictStr.c_str());
      return Action::ERR;
    }
  }
  if (restrictType_ != NONE) {
    cutVal_ = actionArgs.getKeyDouble("cutoff", -1);
    if (cutVal_ <= 0) {
      mprinterr("Error: If 'restrict' is specified, 'cutoff' must be specified and > 0\n");
      return Action::ERR;
    }
  }
    
  property_ = NUMBER;
  if      (actionArgs.hasKey("number") )   property_ = NUMBER;
  else if (actionArgs.hasKey("mass") )     property_ = MASS;
  else if (actionArgs.hasKey("charge") )   property_ = CHARGE;
  else if (actionArgs.hasKey("electron") ) property_ = ELECTRON;

  binType_ = CENTER;
  if      (actionArgs.hasKey("bincenter")) binType_ = CENTER;
  else if (actionArgs.hasKey("binedge")  ) binType_ = EDGE;

  delta_ = actionArgs.getKeyDouble("delta", 0.01);
  if (delta_ <= 0) {
    mprinterr("Error: Delta must be > 0.0\n");
    return Action::ERR;
  }

  // for compatibility with ptraj, ignored because we rely on the atom code to
  // do the right thing, see Atom.{h,cpp}
  if (actionArgs.hasKey("efile"))
    mprintf("Warning: The 'efile' keyword is deprecated.\n");

  // read the rest of the command line as a series of masks
  std::string maskstr;

  unsigned int idx = 1;
  while ( (maskstr = actionArgs.GetMaskNext() ) != emptystring) {
    masks_.push_back( AtomMask(maskstr) );
    if (dsname.empty())
      dsname = init.DSL().GenerateDefaultName("DENSITY");
    MetaData MD(dsname, "avg", idx);
    MD.SetTimeSeries( MetaData::NOT_TS );
    // Hold average density
    DataSet* ads = init.DSL().AddSet( DataSet::DOUBLE, MD );
    if (ads == 0) return Action::ERR;
    ads->SetLegend( NoWhitespace(masks_.back().MaskExpression()) );
    AvSets_.push_back( ads );
    if (outfile != 0) outfile->AddDataSet( ads );
    // Hold SD density
    MD.SetAspect("sd");
    DataSet* sds = init.DSL().AddSet( DataSet::DOUBLE, MD );
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
    // If no masks assume we want total system density.
    if (dsname.empty())
      dsname = actionArgs.GetStringNext();
    density_ = init.DSL().AddSet(DataSet::DOUBLE, dsname, "DENSITY");
    if (density_ == 0) return Action::ERR;
    if (outfile != 0) outfile->AddDataSet( density_ );
    // Hijack delta for storing sum of masses
    delta_ = 0.0;
  } else {
    // Density selected by mask(s) along an axis
    density_ = 0;
    histograms_.resize(masks_.size() );
  }

  mprintf("    DENSITY:");
  if (density_ == 0) {
    const char* binStr[] = {"center", "edge"};
    mprintf(" Determining %s density for %zu masks.\n", PropertyStr_[property_], masks_.size());
    mprintf("\troutine version: %s\n", ROUTINE_VERSION_STRING);
    mprintf("\tDelta is %f\n", delta_);
    mprintf("\tAxis is %s\n", AxisStr_[axis_]);
    mprintf("\tData set name is '%s'\n", dsname.c_str());
    mprintf("\tData set aspect [avg] holds the mean, aspect [sd] holds standard deviation.\n");
    mprintf("\tBin coordinates will be to bin %s.\n", binStr[binType_]);
    if (restrictType_ != NONE) {
      if (restrictType_ == CYLINDER) {
        mprintf("\tOnly atoms within a radius of %f Ang. from the axis will be binned.\n", cutVal_);
        // Square the cutoff
        cutVal_ *= cutVal_;
      } else if (restrictType_ == SQUARE) {
        mprintf("\tOnly atoms within a square of %f Ang. from the axis will be binned.\n", cutVal_);
      } else {
        mprinterr("Internal Error: Restrict type not handled.\n");
        return Action::ERR;
      }
    }
  } else {
    mprintf(" No masks specified, calculating total system density in g/cm^3.\n");
    mprintf("\tData set name is '%s'\n", density_->legend());
  }
  if (outfile != 0)
    mprintf("\tOutput to '%s'\n", outfile->DataFilename().full());

  return Action::OK;
}

// Action_Density::HistSetup()
/** Set up masks and properties for selected atoms. */
Action::RetType Action_Density::HistSetup(ActionSetup& setup) {
  properties_.clear();
  properties_.reserve( masks_.size() );

  int total_nselected = 0;
  for (std::vector<AtomMask>::iterator mask = masks_.begin();
                                       mask != masks_.end();
                                       mask++)
  {
    if (setup.Top().SetupIntegerMask(*mask) ) return Action::ERR;
    mprintf("\t");
    mask->BriefMaskInfo();
    mprintf("\n");
    total_nselected += mask->Nselected();

    std::vector<double> property;

    for (AtomMask::const_iterator idx = mask->begin();
                                  idx != mask->end(); idx++)
    {
      const Atom& atom = setup.Top()[*idx];

      switch (property_) {
        case NUMBER   : property.push_back(1.0); break;
        case MASS     : property.push_back(atom.Mass() ); break;
        case CHARGE   : property.push_back(atom.Charge() ); break;
        case ELECTRON : property.push_back(atom.AtomicNumber() - atom.Charge() ); break;
      }
    }

    properties_.push_back(property);

  }
  if (total_nselected < 1) {
    mprintf("Warning: Nothing selected by masks, skipping.\n");
    return Action::SKIP;
  }

  return Action::OK;  
}

// Action_Density::DensitySetup()
Action::RetType Action_Density::DensitySetup(ActionSetup& setup) {
  // Total system density setup
  if ( !setup.CoordInfo().TrajBox().HasBox() ) { 
    mprintf("Warning: No unit cell information, total density cannot be calculated for '%s'\n",
            setup.Top().c_str());
    return Action::SKIP;
  }
  // delta_ will hold sum of masses
  delta_ = 0.0;
  for (int idx = 0; idx != setup.Top().Natom(); idx++)
    delta_ += setup.Top()[idx].Mass();
  mprintf("\tSum of masses is %g amu\n", delta_);

  return Action::OK;
}

// Action_Density::Setup()
Action::RetType Action_Density::Setup(ActionSetup& setup) {
  if (density_ == 0)
    return HistSetup(setup);
  else
    return DensitySetup(setup);
}

// Action_Density::HistAction()
Action::RetType Action_Density::HistAction(int frameNum, ActionFrame& frm) {
  long int bin = 0;
  // Loop over masks
  for (unsigned int idx = 0; idx != masks_.size(); ++idx)
  {
    AtomMask const& mask = masks_[idx];
    HistType&       hist = histograms_[idx];
    std::map<long int, double> Sum;
    // Loop over mask atoms
    unsigned int midx = 0;
    for (AtomMask::const_iterator atm = mask.begin(); atm != mask.end(); ++atm, midx++)
    {
      const double* XYZ = frm.Frm().XYZ( *atm );
      // Check if we are doing a geometric restriction
      bool binValues = false;
      if (restrictType_ == CYLINDER) {
        //mprintf("DEBUG: CYLINDER RESTRICT: atom %i ordinates %i %i\n", *atm+1, area_coord_[0], area_coord_[1]);
        double x = XYZ[area_coord_[0]];
        double y = XYZ[area_coord_[1]];
        //mprintf("DEBUG:    coords: x= %f  y= %f\n", x, y);
        double d2 = x*x + y*y;
        if (d2 < cutVal_)
          binValues = true;
        //mprintf("DEBUG:    binValues= %i\n", (int)binValues);
      } else if (restrictType_ == SQUARE) {
        if (fabs(XYZ[area_coord_[0]]) < cutVal_ && 
            fabs(XYZ[area_coord_[1]]) < cutVal_)
          binValues = true;
      } else
        binValues = true;
      if (binValues) {
        if (XYZ[axis_] < 0) {
          // Coordinate is negative. Need to subtract off delta so that proper bin
          // is populated (-delta to 0.0).
          bin = (XYZ[axis_] - delta_) / delta_;
        } else {
          // Coordinate is positive.
          bin = XYZ[axis_] / delta_;
        }
        //mprintf("DEBUG: frm=%6i mask=%3u atm=%8i crd=%8.3f bin=%li\n", frameNum+1, idx, *atm+1, XYZ[axis_], bin);
        Sum[bin] += properties_[idx][midx];
      }
    } // END loop over atoms in mask
    // Accumulate sums
    if (!Sum.empty())
      hist.accumulate( Sum );
  }

  // Accumulate area
  Box const& box = frm.Frm().BoxCrd();
  area_.accumulate(box.Param((Box::ParamType)area_coord_[0]) *
                   box.Param((Box::ParamType)area_coord_[1]));

  return Action::OK;
}

/** Convert units from amu/Ang^3 to g/cm^3 */
const double Action_Density::AMU_ANG_TO_G_CM3 = Constants::NA * 1E-24;

// Action_Density::DensityAction()
Action::RetType Action_Density::DensityAction(int frameNum, ActionFrame& frm) {
  // Total mass is in delta_
  double density = delta_ / (frm.Frm().BoxCrd().CellVolume() * AMU_ANG_TO_G_CM3);
  density_->Add(frameNum, &density);
  return Action::OK;
}

// Action_Density::DoAction()
Action::RetType Action_Density::DoAction(int frameNum, ActionFrame& frm) {
  if (density_ == 0)
    return HistAction(frameNum, frm);
  else
    return DensityAction(frameNum, frm);
}

#ifdef MPI
/** Combine histogram data. */
int Action_Density::SyncAction() {
  if (trajComm_.Size() < 2) return 0;
  std::vector<double> dbuffer;
  std::vector<long int> ibuffer;
  unsigned long rank_size;
  // Should always be same number of histograms on each rank
  if (trajComm_.Master()) {
    // Master.
    for (int rank = 1; rank < trajComm_.Size(); rank++)
    {
      for (HistArray::iterator hist = histograms_.begin(); hist != histograms_.end(); ++hist)
      {
        // 1. Receive number of bins for hist from rank
        //rprintf("DEBUG: Receiving bins from rank %i\n", rank);
        trajComm_.SendMaster(&rank_size, 1, rank, MPI_UNSIGNED_LONG);
        //rprintf("DEBUG: %lu bins on rank %i\n", rank_size, rank);
        // 2. Recieve histogram Keys and values from rank
        // Double values: n, mean[], m2[]
        dbuffer.resize( 1 + (2*rank_size) );
        ibuffer.resize( rank_size );
        trajComm_.SendMaster(&ibuffer[0], rank_size,       rank, MPI_LONG);
        trajComm_.SendMaster(&dbuffer[0], 1+(2*rank_size), rank, MPI_DOUBLE);
        // Combine histograms
        hist->Combine( HistType(ibuffer, dbuffer) );
      } // END loop over master histograms
    } // END loop over ranks
  } else {
    // Not master. Send histogram data to master.
    for (HistArray::const_iterator hist = histograms_.begin(); hist != histograms_.end(); ++hist)
    {
      // 1. Send number of bins for hist on this rank to master
      rank_size = (unsigned long)hist->nBins();
      trajComm_.SendMaster(&rank_size, 1, trajComm_.Rank(), MPI_UNSIGNED_LONG);
      // 2. Send histogram Keys and values to master.
      ibuffer.clear();
      dbuffer.clear();
      // Double values: n, mean[], m2[]
      dbuffer.reserve( 1 + (2*hist->nBins()) );
      ibuffer.reserve( hist->nBins() );
      dbuffer.push_back( hist->nData() );
      for (HistType::const_iterator it = hist->mean_begin(); it != hist->mean_end(); ++it) {
        ibuffer.push_back( it->first );
        dbuffer.push_back( it->second );
      }
      for (HistType::const_iterator it = hist->variance_begin(); it != hist->variance_end(); ++it)
        dbuffer.push_back( it->second );
      trajComm_.SendMaster(&ibuffer[0], hist->nBins(),       trajComm_.Rank(), MPI_LONG);
      trajComm_.SendMaster(&dbuffer[0], 1+(2*hist->nBins()), trajComm_.Rank(), MPI_DOUBLE);
    } // END loop over rank histograms
  }
  return 0;
}
#endif /* MPI */

// Action_Density::PrintHist()
/** Do histogram normalization. */
void Action_Density::PrintHist()
{
  // Determine if area scaling should occur
  const unsigned int SMALL = 1.0;

  double avgArea = area_.mean();
  double sdArea  = sqrt(area_.variance());
  bool scale_area = (property_ == ELECTRON && avgArea > SMALL);

  mprintf("    DENSITY: The average box area in %c/%c is %.2f Angstrom (sd = %.2f).\n",
          area_coord_[0] + 88, area_coord_[1] + 88, avgArea, sdArea);

  if (scale_area)
    mprintf("The electron density will be scaled by this area.\n");

  // Loop over all histograms. Find lowest and highest bin idx out of all
  // histograms. All output histograms will have the same dimensions.
  long int lowest_idx = 0;
  long int highest_idx = 0;
  for (unsigned int idx = 0; idx != histograms_.size(); idx++) 
  {
    HistType const& hist = histograms_[idx];
    if (hist.empty()) {
      mprintf("Warning: Histogram for '%s' is empty; skipping.\n", masks_[idx].MaskString());
      continue;
    }
    // Find lowest bin
    if (idx == 0) {
      lowest_idx  = hist.lowestKey();
      highest_idx = hist.highestKey();
    } else {
      lowest_idx  = std::min(lowest_idx,  hist.lowestKey());
      highest_idx = std::max(highest_idx, hist.highestKey());
    }
  }
  // If using center of bins, put blank bins at either end.
  double xshift = 0.0;
  if (binType_ == CENTER) {
    lowest_idx--;
    highest_idx++;
    xshift = delta_ / 2;
  }
  // Set up common output dimensions
  long int Nbins = (highest_idx - lowest_idx + 1);
  long int offset = -lowest_idx;
  double Xmin = (delta_ * lowest_idx) + xshift;
  //mprintf("DEBUG: Lowest idx= %li  Xmin= %g  Highest idx= %li  Bins= %li\n",
  //        lowest_idx, Xmin, highest_idx, Nbins);
  Dimension Xdim(Xmin, delta_, AxisStr_[axis_]);
  // Loop over all histograms. Normalize and populate output sets.
  for (unsigned int idx = 0; idx != histograms_.size(); idx++) 
  {
    HistType const& hist = histograms_[idx];
    if (hist.empty()) {
      continue;
    }
    // Calculate normalization
    double fac   = delta_;
    double sdfac = 1.0;
    if (scale_area) {
      fac *= avgArea;
      sdfac = 1.0 / avgArea;
    }
    fac = 1.0 / fac;
    // Populate output data sets
    DataSet_1D& out_av = static_cast<DataSet_1D&>( *(AvSets_[idx]) );
    DataSet_1D& out_sd = static_cast<DataSet_1D&>( *(SdSets_[idx]) );
    out_av.Allocate(DataSet::SizeArray(1, Nbins));
    out_sd.Allocate(DataSet::SizeArray(1, Nbins));
    out_av.SetDim(Dimension::X, Xdim);
    out_sd.SetDim(Dimension::X, Xdim);
    // Blank bin for bin center
    if (binType_ == CENTER) {
      double dzero = 0.0;
      out_av.Add(0, &dzero);
      out_sd.Add(0, &dzero);
    }
    // Loop over populated bins
    HistType::const_iterator var = hist.variance_begin();
    for (HistType::const_iterator mean = hist.mean_begin();
                                  mean != hist.mean_end(); ++mean, ++var)
    {
      long int frm = mean->first + offset;
      double density  = mean->second * fac;
      out_av.Add(frm, &density);
      double variance;
      if (hist.nData() < 2)
        variance = 0;
      else
        variance = var->second / (hist.nData() - 1);
      if (variance > 0) {
        variance = sqrt(variance);
        variance *= sdfac;
      }
      out_sd.Add(frm, &variance);
    }
    // Blank bin for bin center
    if (binType_ == CENTER) {
      double dzero = 0.0;
      out_av.Add(Nbins-1, &dzero);
      out_sd.Add(Nbins-1, &dzero);
    }
  } // END loop over all histograms
}

// Action_Density::PrintDensity()
void Action_Density::PrintDensity() {
  if (density_->Size() > 0) {
    double stdev;
    double avg = ((DataSet_1D*)density_)->Avg(stdev);
    mprintf("    DENSITY: Avg= %g  Stdev= %g (%zu elements), g/cm^3\n",
            avg, stdev, density_->Size());
  }
}

// Action_Density::Print()
void Action_Density::Print() {
  if (density_ == 0)
    PrintHist();
  else
    PrintDensity();
}
