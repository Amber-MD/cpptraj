#include "Analysis_Hist.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // doubleToString
#include "DS_Math.h" // Min, Max
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
int Analysis_Hist::setupDimension(ArgList &arglist, DataSet_1D const& dset) {
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
  if (!minArg) {
    if (!minArgSet_) 
      dim.SetMin( DS_Math::Min(dset) );
    else
      dim.SetMin( default_dim_.Min() );
  }
  // If no max arg and no default max arg, get max from dataset
  if (!maxArg) {
    if (!maxArgSet_)
      dim.SetMax( DS_Math::Max(dset) );
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
  int offset = 1;
  for (std::vector<Dimension>::reverse_iterator rd = dimensions_.rbegin();
                                                rd != dimensions_.rend(); ++rd)
  {
    if (debug_>0) mprintf("\tHistogram: %s offset is %i\n",(*rd).Label().c_str(), offset);
    (*rd).SetOffset( offset );
    offset *= (*rd).Bins();
  }
  // offset should now be equal to the total number of bins across all dimensions
  if (debug_>0) mprintf("\tHistogram: Total Bins = %i\n",offset);

  return 0;
}

// Analysis_Hist::Setup()
/** Set up histogram with specified data sets. */
Analysis::RetType Analysis_Hist::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  debug_ = debugIn;

  // Keywords
  std::string histname = analyzeArgs.GetStringKey("name");
  std::string outfilename = analyzeArgs.GetStringKey("out");
  if (outfilename.empty()) {
    mprintf("Error: Hist: No output filename specified.\n");
    return Analysis::ERR;
  }
  Temp_ = analyzeArgs.getKeyDouble("free",-1.0);
  if (Temp_!=-1.0) calcFreeE_ = true;
  //gnuplot_ = analyzeArgs.hasKey("gnu");
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
      mprinterr("Internal Error: Currently histogram beyond 3 dimensions disabled.\n");
      return Analysis::ERR;
  }
  // Set up output data file
  outfile_ = DFLin->AddDataFile(outfilename, analyzeArgs);
  if (outfile_==0) return Analysis::ERR;
  
  mprintf("\tHist: %s: Set up for %zu dimensions using the following datasets:\n", 
          outfilename.c_str(), N_dimensions_);
  mprintf("\t      [ ");
  for (std::vector<DataSet_1D*>::iterator ds=histdata_.begin(); ds!=histdata_.end(); ++ds)
    mprintf("%s ",(*ds)->Legend().c_str());
  mprintf("]\n");
  if (calcFreeE_)
    mprintf("\t      Free energy will be calculated from bin populations at %lf K.\n",Temp_);
  /*if (circular_ || gnuplot_) {
    mprintf("\tWarning: gnuplot and/or circular specified; advanced grace/gnuplot\n");
    mprintf("\t         formatting disabled.\n");
    if (circular_)
      mprintf("\t      circular: Output coordinates will be wrapped.\n");
    if (gnuplot_)
      mprintf("\t      gnuplot: Output will be in gnuplot-readable format.\n");
  }*/
  if (normalize_)
    mprintf("\t      norm: Bins will be normalized to 1.0.\n");

  return Analysis::OK;
}

// Analysis_Hist::Analyze()
Analysis::RetType Analysis_Hist::Analyze() {
  // Set up dimensions
  // Size of histdata and dimensionArgs should be the same
  for (unsigned int hd = 0; hd < N_dimensions_; hd++) {
    if ( setupDimension(dimensionArgs_[hd], *(histdata_[hd])) ) 
      return Analysis::ERR;
  }

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

  // Set up appropriate data set type for binning
  DataSet_double* H1D = 0;
  DataSet_MatrixDbl* H2D = 0;
  DataSet_GridFlt* H3D = 0;
  switch (N_dimensions_) {
    case 1: H1D = (DataSet_double*)hist_; break;
    case 2: H2D = (DataSet_MatrixDbl*)hist_; break;
    case 3: H3D = (DataSet_GridFlt*)hist_; break;
  }
  // Bin data
  for (size_t n = 0; n < Ndata; n++) {
    int index = 0;
    std::vector<Dimension>::iterator dim = dimensions_.begin();
    for (std::vector<DataSet_1D*>::iterator ds = histdata_.begin();
                                            ds != histdata_.end(); ++ds)
    {
      double dval = (*ds)->Dval( n );
      // Check if data is out of bounds.
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
      switch (N_dimensions_) {
        case 1: (*H1D)[index]++; break;
        case 2: (*H2D)[index]++; break;
        case 3: (*H3D)[index]++; break;
      }
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
  if (calcFreeE_) hist_->CalcFreeE(Temp_, -1);

  // Normalize if requested
  if (normalize_) hist_->Normalize();

  if (outfile_ == 0) {
    // Use Histogram built-in output
    hist_->PrintBins(outfilename_.c_str(), circular_, gnuplot_);
  } else {
    outfile_->AddSet( hist_ );
    // Using DataFileList framework, set-up labels etc.
    if (hist_->NumDimension() == 1) {
      outfile_->ProcessArgs("xlabel " + (*hist_)[0].Label() + " ylabel Count" +
                            " xmin "  + doubleToString( (*hist_)[0].Min()  ) +
                            " xstep " + doubleToString( (*hist_)[0].Step() )  );
    } else { // two dimensions
      outfile_->ProcessArgs("xlabel "  + (*hist_)[0].Label() +
                            " ylabel " + (*hist_)[1].Label() +
                            " xmin "   + doubleToString( (*hist_)[0].Min()  ) +
                            " xstep " + doubleToString( (*hist_)[0].Step() ) +
                            " ymin "  + doubleToString( (*hist_)[1].Min()  ) +
                            " ystep " + doubleToString( (*hist_)[1].Step() ) +
                            "noxcol usemap nolabels");
    }
  }

  //hist.PrintBins(false,false);
  return Analysis::OK;
}
