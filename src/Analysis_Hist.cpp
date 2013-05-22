#include <cmath> // log
#include "Analysis_Hist.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // doubleToString
// DataSet types used by Analysis_Hist
#include "DataSet_double.h"
#include "DataSet_MatrixDbl.h"
#include "DataSet_GridFlt.h"

// CONSTRUCTOR
Analysis_Hist::Analysis_Hist() :
  outfile_(0), 
  hist_(0),
  debug_(0), 
  calcFreeE_(false),
  Temp_(-1.0),
  normalize_(false),
//  gnuplot_(false),
  circular_(false),
  N_dimensions_(0),
  minArgSet_(false),
  maxArgSet_(false)
{}

void Analysis_Hist::Help() {
  mprintf("\t<dataset_name>[,min,max,step,bins] ...\n");
  mprintf("\t[free <temperature>] [norm] [gnu] [circular] out <filename>\n");
  mprintf("\t[min <min>] [max <max>] [step <step>] [bins <bins>]\n");
  mprintf("\tHistogram the given data set(s)\n");
}

// Analysis_Hist::CheckDimension()
/** Given an argument with format, DataSet_Name[,min,max,step,bins], check
  * that DataSet_Name exists and is valid. Add the argument to 
  * dimensionArgs and the corresponding dataset to histdata.
  */
int Analysis_Hist::CheckDimension(std::string const& input, DataSetList *datasetlist) {
  ArgList arglist;
  // Separate input string by ','
  arglist.SetList(input, ",");
  if (arglist.Nargs()<1) {
    mprintf("Warning: Hist::CheckDimension: No arguments found in input: %s\n",input.c_str());
    return 1;
  }

  // First argument should specify dataset name
  if (debug_>0) mprintf("\tHist: Setting up histogram dimension using dataset %s\n",
                       arglist.Command());
  DataSet* dset = datasetlist->GetDataSet( arglist[0] );
  if (dset == 0) {
    mprintf("\t      Dataset %s not found.\n",arglist.Command());
    return 1;
  }

  // For now only 1D data sets can be histogrammed
  if (dset->Ndim() != 1) {
    mprinterr("Error: Hist: dataset %s has %u dimensions.\n",
              dset->Legend().c_str(), dset->Ndim());
    mprinterr("Error: Hist: Currently only 1D data sets can be histogrammed.\n");
    return 1;
  }

  // Check that dataset is not string
  if (dset->Type()==DataSet::STRING) {
    mprintf("Error: Hist: Cannot histogram dataset %s, type STRING.\n", 
            dset->Legend().c_str());
    return 1;
  }

  dimensionArgs_.push_back( arglist );
  histdata_.push_back( (DataSet_1D*)dset );
  return 0;
}

// Analysis_Hist::setupDimension()
/** Given an ArgList containing name,[min,max,step,bins,col,N], set up a 
  * coordinate with that name and parameters min, max, step, bins.
  * If '*' or not specified, a default value will be set.
  * \return 1 if error occurs, 0 otherwise.
  */
int Analysis_Hist::setupDimension(ArgList &arglist, DataSet_1D const& dset, size_t& offset) {
  Dimension dim;
  bool minArg = false;
  bool maxArg = false;

  if (debug_>1)
    arglist.PrintList();

  // Set up dimension name
  // NOTE: arglist[0] should be same as dset name from CheckDimension 
  dim.SetLabel( arglist[0] );

  // Cycle through coordinate arguments. Any argument left blank will be 
  // assigned a default value later.
  for (int i=1; i<arglist.Nargs(); i++) {
    if (debug_>1) mprintf("    DEBUG: setupCoord: Token %i (%s)\n",i,arglist[i].c_str());
    // Default explicitly requested
    if (arglist[i] == "*") continue;
    switch (i) {
      case 1 : dim.SetMin(  convertToDouble( arglist[i]) ); minArg=true; break;
      case 2 : dim.SetMax(  convertToDouble( arglist[i]) ); maxArg=true; break;
      case 3 : dim.SetStep( convertToDouble( arglist[i]) ); break;
      case 4 : dim.SetBins( convertToInteger(arglist[i]) ); break;
    }
  }

  // If no min arg and no default min arg, get min from dataset
  // TODO: Use Min/MaxIsSet
  if (!minArg) {
    if (!minArgSet_) 
      dim.SetMin( dset.Min() );
    else
      dim.SetMin( default_dim_.Min() );
  }
  // If no max arg and no default max arg, get max from dataset
  if (!maxArg) {
    if (!maxArgSet_)
      dim.SetMax( dset.Max() );
    else
      dim.SetMax( default_dim_.Max() );
  }
  // Check that min < max
  if (dim.Min() >= dim.Max()) {
    mprinterr("Error: Hist: Dimension %s: min (%lf) must be less than max (%lf).\n",
              dim.Label().c_str(), dim.Min(), dim.Max());
    return 1;
  }

  // If bins/step not specified, use default
  if (dim.Bins()==-1)
    dim.SetBins( default_dim_.Bins() );
  if (dim.Step()==-1)
    dim.SetStep( default_dim_.Step() );

  // Attempt to set up bins or step.
  if (dim.CalcBinsOrStep()!=0) return 1;
 
  dim.PrintDim();
  dimensions_.push_back( dim );

  // Recalculate offsets for all dimensions starting at farthest coord. This
  // follows column major ordering.
  offset = 1UL;
  for (std::vector<Dimension>::reverse_iterator rd = dimensions_.rbegin();
                                                rd != dimensions_.rend(); ++rd)
  {
    if (debug_>0) mprintf("\tHistogram: %s offset is %zu\n",(*rd).Label().c_str(), offset);
    (*rd).SetOffset( offset );
    offset *= (*rd).Bins();
  }
  // offset should now be equal to the total number of bins across all dimensions
  if (debug_>0) mprintf("\tHistogram: Total Bins = %zu\n",offset);

  return 0;
}

// Analysis_Hist::Setup()
/** Set up histogram with specified data sets. */
Analysis::RetType Analysis_Hist::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  //debug_ = debugIn; // DEBUG
  debug_ = 1; // DEBUG - enable debug
  mprintf("Warning: HIST DEBUG ACTIVE\n"); // DEBUG

  // Keywords
  std::string histname = analyzeArgs.GetStringKey("name");
  outfilename_ = analyzeArgs.GetStringKey("out");
  if (outfilename_.empty()) {
    mprintf("Error: Hist: No output filename specified.\n");
    return Analysis::ERR;
  }
  Temp_ = analyzeArgs.getKeyDouble("free",-1.0);
  if (Temp_!=-1.0) calcFreeE_ = true;
  gnuplot_ = analyzeArgs.hasKey("gnu");
  normalize_ = analyzeArgs.hasKey("norm");
  circular_ = analyzeArgs.hasKey("circular");
  if ( analyzeArgs.Contains("min") ) {
    default_dim_.SetMin( analyzeArgs.getKeyDouble("min",0.0) );
    minArgSet_ = true;
  }
  if ( analyzeArgs.Contains("max") ) {
    default_dim_.SetMax( analyzeArgs.getKeyDouble("max",0.0) );
    maxArgSet_ = true;
  }
  default_dim_.SetStep( analyzeArgs.getKeyDouble("step",-1.0) );
  default_dim_.SetBins( analyzeArgs.getKeyInt("bins",-1) );

  // Treat all remaining arguments as dataset names. Do not set up dimensions
  // yet since the data sets may not be fully populated.
  ArgList dsetNames = analyzeArgs.RemainingArgs();
  for ( ArgList::const_iterator setname = dsetNames.begin(); 
                                setname != dsetNames.end(); ++setname)
  { 
    if (CheckDimension( *setname, datasetlist )) return Analysis::ERR;
  }
  // histdata contains the DataSets to be histogrammed
  if (histdata_.empty()) {
    mprinterr("Error: Hist: No datasets specified.\n");
    return Analysis::ERR;
  }
  // Total # of dimensions for the histogram is the number of sets to be binned.
  N_dimensions_ = histdata_.size();
  switch ( N_dimensions_ ) {
    case 1: hist_ = datasetlist->AddSet( DataSet::DOUBLE,     histname, "Hist"); break;
    case 2: hist_ = datasetlist->AddSet( DataSet::MATRIX_DBL, histname, "Hist"); break;
    // TODO: GRID_DBL
    case 3: hist_ = datasetlist->AddSet( DataSet::GRID_FLT,   histname, "Hist"); break;
    default: // FIXME: GET N DIMENSION CASE!
      mprintf("Warning: Histogram dimension > 3. DataSet/DataFile output not supported.\n");
      mprintf("Warning: Using internal routine for output.\n");
  }
  // Set up output data file
  outfile_ = 0;
  if (N_dimensions_ < 3) {
    if (hist_ == 0) {
      mprinterr("Error: Could not set up histogram data set.\n");
      return Analysis::ERR;
    }
    outfile_ = DFLin->AddDataFile(outfilename_, analyzeArgs);
    if (outfile_==0) return Analysis::ERR;
    outfile_->AddSet( hist_ );
  }
  
  mprintf("\tHist: %s: Set up for %zu dimensions using the following datasets:\n", 
          outfilename_.c_str(), N_dimensions_);
  mprintf("\t      [ ");
  for (std::vector<DataSet_1D*>::iterator ds=histdata_.begin(); ds!=histdata_.end(); ++ds)
    mprintf("%s ",(*ds)->Legend().c_str());
  mprintf("]\n");
  if (calcFreeE_)
    mprintf("\t      Free energy will be calculated from bin populations at %lf K.\n",Temp_);
  //if (circular_ || gnuplot_) {
  //  mprintf("\tWarning: gnuplot and/or circular specified; advanced grace/gnuplot\n");
  //  mprintf("\t         formatting disabled.\n");*/
    if (circular_)
      mprintf("\t      circular: Output coordinates will be wrapped.\n");
    if (gnuplot_ && outfile_ == 0)
      mprintf("\t      gnuplot: Output will be in gnuplot-readable format.\n");
  //}
  if (normalize_)
    mprintf("\t      norm: Bins will be normalized to 1.0.\n");

  return Analysis::OK;
}

// Analysis_Hist::Analyze()
Analysis::RetType Analysis_Hist::Analyze() {
  // Set up dimensions
  // Size of histdata and dimensionArgs should be the same
  size_t total_bins = 0UL;
  for (unsigned int hd = 0; hd < N_dimensions_; hd++) {
    if ( setupDimension(dimensionArgs_[hd], *(histdata_[hd]), total_bins) ) 
      return Analysis::ERR;
  }
  // dimensionArgs no longer needed
  dimensionArgs_.clear();

  // Check that the number of data points in each dimension are equal
  std::vector<DataSet_1D*>::iterator ds = histdata_.begin();
  size_t Ndata = (*ds)->Size();
  ++ds;
  for (; ds != histdata_.end(); ++ds)
  {
    //mprintf("DEBUG: DS %s size %i\n",histdata[hd]->Name(),histdata[hd]->Xmax()+1);
    if (Ndata != (*ds)->Size()) {
      mprinterr("Error: Hist: Dataset %s has inconsistent # data points (%u), expected %u.\n",
                (*ds)->Legend().c_str(), (*ds)->Size(), Ndata);
      return Analysis::ERR;
    }
  }
  mprintf("\tHist: %u data points in each dimension.\n", Ndata);

  // Allocate bins
  mprintf("\tHist: Allocating histogram, total bins = %zu\n", total_bins);
  Bins_.resize( total_bins, 0.0 );

  // Bin data
  for (size_t n = 0; n < Ndata; n++) {
    int index = 0; // TODO: Make long?
    std::vector<Dimension>::iterator dim = dimensions_.begin();
    for (std::vector<DataSet_1D*>::iterator ds = histdata_.begin();
                                            ds != histdata_.end(); ++ds, ++dim)
    {
      double dval = (*ds)->Dval( n );
      // Check if data is out of bounds for this dimension.
      if (dval > (*dim).Max() || dval < (*dim).Min()) {
        index = -1;
        break;
      }
      // Calculate index for this particular dimension (idx)
      int idx = (int)((dval - (*dim).Min()) / (*dim).Step());
      if (debug_>1) mprintf(" [%s:%f (%i)],",(*dim).Label().c_str(), dval, idx);
      // Calculate overall index in Bins, offset has already been calcd.
      index += (idx * (*dim).Offset());
    }
    // If index was successfully calculated, populate bin
    if (index > -1 && index < (int)hist_->Size()) {
      if (debug_ > 1) mprintf(" |index=%i",index);
      Bins_[index]++;
    } else {
      mprintf("\tWarning: Frame %u Coordinates out of bounds (%i)\n", n, index);
    }
    if (debug_>1) mprintf("}\n");
  }
/*
  std::vector<double> coord( hist_->NumDimension() );
  for (int n=0; n < Ndata; n++) {
    std::vector<double>::iterator coord_it = coord.begin();
    for (std::vector<DataSet*>::iterator ds = histdata_.begin(); ds != histdata_.end(); ++ds) {
      *coord_it = (*ds)->Dval( n );
      ++coord_it;
    }
    hist_->BinData( coord );
  }
*/
  // Calc free energy if requested
  if (calcFreeE_) CalcFreeE();

  // Normalize if requested
  if (normalize_) Normalize();

  // Using DataFileList framework, set-up labels etc.
  if (N_dimensions_ == 1) {
    DataSet_double& dds = static_cast<DataSet_double&>( *hist_ );
    dds.Allocate1D( dimensions_[0].Bins() );
    std::copy( Bins_.begin(), Bins_.end(), dds.begin() );
    outfile_->ProcessArgs("xlabel " + dimensions_[0].Label() + " ylabel Count" +
                          " xmin "  + doubleToString( dimensions_[0].Min()  ) +
                          " xstep " + doubleToString( dimensions_[0].Step() )  );
  } else if (N_dimensions_ == 2) {
    DataSet_MatrixDbl& mds = static_cast<DataSet_MatrixDbl&>( *hist_ );
    mds.Allocate2D( dimensions_[0].Bins(), dimensions_[1].Bins() );
    std::copy( Bins_.begin(), Bins_.end(), mds.begin() );
    outfile_->ProcessArgs("xlabel "  + dimensions_[0].Label() +
                          " ylabel " + dimensions_[1].Label() +
                          " xmin "   + doubleToString( dimensions_[0].Min()  ) +
                          " xstep "  + doubleToString( dimensions_[0].Step() ) +
                          " ymin "   + doubleToString( dimensions_[1].Min()  ) +
                          " ystep "  + doubleToString( dimensions_[1].Step() ) +
                          "noxcol usemap nolabels");
  } else if (N_dimensions_ == 3) {
    DataSet_GridFlt& gds = static_cast<DataSet_GridFlt&>( *hist_ );
    gds.Allocate3D( dimensions_[0].Bins(), dimensions_[1].Bins(), dimensions_[2].Bins() );
    std::copy( Bins_.begin(), Bins_.end(), gds.begin() );
    outfile_->ProcessArgs("xlabel "  + dimensions_[0].Label() +
                          " ylabel " + dimensions_[1].Label() +
                          " zlabel " + dimensions_[2].Label() + 
                          " xmin "   + doubleToString( dimensions_[0].Min()  ) +
                          " xstep "  + doubleToString( dimensions_[0].Step() ) +
                          " ymin "   + doubleToString( dimensions_[1].Min()  ) +
                          " ystep "  + doubleToString( dimensions_[1].Step() ) +
                          " zmin "   + doubleToString( dimensions_[2].Min()  ) +
                          " zstep "  + doubleToString( dimensions_[2].Step() ) +
                          "noxcol usemap nolabels");
  } else {
    // Use Histogram built-in output
    PrintBins();
  }

  //hist.PrintBins(false,false);
  return Analysis::OK;
}

/** Calculate free energy based on bin populations.  */
int Analysis_Hist::CalcFreeE() {
  int refbin = -1; // TODO: re-enable
    mprintf("\tHistogram: Calculating free E at %f K.\n",Temp_);
  // TODO: Make Kb a constant
  double KT = (-.001985886596 * Temp_);

  // Find most populated bin for G=0
  std::vector<double>::iterator bin = Bins_.begin();
  double binmax = *bin;
  ++bin;
  for (; bin != Bins_.end(); ++bin)
    if (*bin > binmax)
      binmax = *bin;
  mprintf("\t           Bins max is %.0f\n",binmax);
  if (binmax==0) {
    mprinterr("Histogram: Cannot calc free E, no bins populated!\n");
    return 1;
  }

  // If requested, set up reference bin other than max
  if (refbin>-1) {
    if (Bins_[refbin] > 0) {
      binmax = Bins_[refbin];
      mprintf("\t           Calculating free E w.r.t bin %i, population %f\n",refbin,binmax);
    } else
      mprintf("Warning: Reference bin %i has no population. Using %f\n",refbin,binmax);
  }

  // Set artificial ceiling for bins with 0 population. Make it equivalent
  // to a bin population of 0.5. 
  double temp = 0.5 / binmax;
  double ceiling = log(temp) * KT;
  //ceiling*=1.10;
  mprintf("\t           Artificial ceiling (bin pop = 0.5) is %f kcal/mol.\n",ceiling);

  // Calculate free E based on populations
  for (bin = Bins_.begin(); bin != Bins_.end(); ++bin) {
    temp = *bin;               // Store Bin population in temp
    if (temp>0) {
      temp /= binmax;          // bin/max
      *bin = log(temp) * KT;   // -R*T*ln(bin/max)
    } else
      *bin = ceiling;          // Artificial ceiling for 0 pop bins.
  }

  return 0;
}

/** Normalize bins so that sum over all bins is 1.0 */
int Analysis_Hist::Normalize() {
  double sum = 0.0;
  mprintf("\tHistogram: Normalizing bin populations to 1.0\n");
  for (std::vector<double>::iterator bin = Bins_.begin(); bin != Bins_.end(); ++bin)
    sum += *bin;
  mprintf("\t           Sum over all bins is %lf\n",sum);
  if (sum == 0.0) {
    mprinterr("Error: Histogram::Normalize: Sum over bin populations is 0.0\n");
    return 1;
  }
  for (std::vector<double>::iterator bin = Bins_.begin(); bin != Bins_.end(); ++bin)
    *bin /= sum;
  return 0;
}

// -----------------------------------------------------------------------------
/** Set current bin to 0 and initialize indices. If isCircularIn is true the
  * bin indices will wrap in each dimension.
  */
static inline std::vector<int> BinStart(int Ndim, bool circularIn) {
  std::vector<int> BinIndices;
  if (circularIn)
    BinIndices.assign( Ndim, -1 );
  else
    BinIndices.assign( Ndim,  0 );
  //mprintf("BinStart: Dimension %zu, circular is %i\n",BinIndices_.size(), isCircular_);
  //mprintf("D1=%i  D2=%i\n",dimensions_[0].Bins(), dimensions_[1].Bins());
  return BinIndices;
}

/** Calculate index based on BinIndices, converting wrapped indices to 
  * actual indices.
  */
int Analysis_Hist::BinIndicesToIndex(std::vector<int> const& BinIndices) {
  int idx;
  int index = 0;
  std::vector<int>::const_iterator count = BinIndices.begin();
  for (std::vector<Dimension>::iterator dim = dimensions_.begin();
                                        dim != dimensions_.end(); ++dim, ++count)
  {
    //mprinterr(" %i",count[n]);
    if (*count == -1)
      idx = (*dim).Bins() - 1;
    else if (*count == (*dim).Bins())
      idx = 0;
    else
      idx = *count;
    index += (idx * (*dim).Offset());
  }
  //mprinterr(" = %i\n",index);
  return index;
}

// Histogram::IncrementBinIndices()
/** \return true if there are more bins to process.
  * \return false if there are no more bins.
  */
bool Analysis_Hist::IncrementBinIndices(std::vector<int>& BinIndices, 
                                        int isCircular, bool& hasCycled)
{
  //mprintf("DEBUG0\t\t\tCoord0=%i Coord1=%i\n",BinIndices_[0],BinIndices_[1]);
  // Increment highest order coord.
  std::vector<int>::reverse_iterator rcount = BinIndices.rbegin();
  ++(*rcount);
  std::vector<Dimension>::reverse_iterator rdim = dimensions_.rbegin();
  // Check if highest order coord has Cycled.
  if (*rcount == (*rdim).Bins() + isCircular)
    hasCycled = true;
  else
    hasCycled = false;
  // Increment other coords if necessary
  for (; rdim != dimensions_.rend() - 1; ++rdim) {
    if (*rcount == (*rdim).Bins()+isCircular) {
      (*rcount) = 0 - isCircular;
      ++rcount;
      ++(*rcount);
    }
  }
  //mprintf("DEBUG1\t\t\tCoord0=%i Coord1=%i\n",BinIndices_[0],BinIndices_[1]);
  // If the lowest order coord is at lowest dimensions #bins, done.
  if (BinIndices[0] == dimensions_[0].Bins() + isCircular)
    return false;

  return true;
}

/** This routine will print out the contents of the Bin array in column-major 
  * order. A counter is used to calculate the appropriate coordinate indices 
  * for the array. The count array works like this:
  *   Dim0index Dim1index ... DimNindex
  * DimNindex is considered the highest order dimension and is always 
  * incremented. When DimNindex reaches the number of bins for that Dim
  * it cycles, and the next highest order dim is checked, down to the lowest
  * order dimension.
  * If circular specified, wrapping will occur so that data from indices -1 
  * and N+1 (corresponding to N and 0) are printed out as well.
  * --- CURRENTLY NOT IMPLEMENTED ---
  * If binType is 0, the Bins array will be printed. If binType is 1, the
  * landscape will be printed.
  * if SD is not null, the standard deviation array will also be printed.
  */
void Analysis_Hist::PrintBins() {
  int isCircular = 0;
  bool hasCycled = false;
  CpptrajFile outfile;

  if (outfile.SetupWrite(outfilename_, debug_)) return;
  if (outfile.OpenFile()) return;

  mprintf("\tHistogram: Writing standard histogram file %s\n",outfilename_.c_str());

  std::vector<int> BinIndices = BinStart( dimensions_.size(), circular_ );
  if (circular_) isCircular = 1;

  if (gnuplot_) {
    if (dimensions_.size() == 2)
      outfile.Printf("set pm3d map\nsplot \"-\" with pm3d title \"%s\"\n",
                     outfilename_.c_str());
    else if (dimensions_.size() == 1)
      outfile.Printf("plot \"-\"\n",outfilename_.c_str());
  }

  if (debug_>0) {
    if (circular_)
      mprintf("\t\tPrinting %zu bins in circular fashion.\n", Bins_.size());
    else
      mprintf("\t\tPrinting %zu bins.\n",Bins_.size());
  }

  bool loop=true;
  while (loop) {
    int index = BinIndicesToIndex(BinIndices);
    // If we dont care about zero bins or bin pop > 0, output
    for (unsigned int i=0; i < dimensions_.size(); ++i)
      outfile.Printf("%lf ",
                      ((double)BinIndices[i]*dimensions_[i].Step()) + dimensions_[i].Min() );
    outfile.Printf("%lf\n",Bins_[index]);

    loop = IncrementBinIndices(BinIndices, isCircular, hasCycled);
    // If gnuplot, print extra space when highest order coord cycles
    if (gnuplot_ && hasCycled)
        outfile.Printf("\n");
  }
  if (gnuplot_ && dimensions_.size() < 3)
    outfile.Printf("end\npause -1\n");
  outfile.CloseFile();
}

