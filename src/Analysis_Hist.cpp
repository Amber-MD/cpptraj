#include "Analysis_Hist.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // doubleToString
#include "DS_Math.h" // Min, Max
// Analysis_Hist

// CONSTRUCTOR
Analysis_Hist::Analysis_Hist() :
  outfile_(0), 
  hist_(0),
  debug_(0), 
  calcFreeE_(false),
  Temp_(-1.0),
  normalize_(false),
  gnuplot_(false),
  circular_(false),
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
  DataSet *dset = datasetlist->GetDataSet( arglist[0] );
  if (dset == 0) {
    mprintf("\t      Dataset %s not found.\n",arglist.Command());
    return 1;
  }

  // Check that dataset is not string
  if (dset->Type()==DataSet::STRING) {
    mprintf("Error: Hist: Cannot histogram dataset %s, type STRING.\n", 
            dset->Legend().c_str());
    return 1;
  }

  dimensionArgs_.push_back( arglist );
  histdata_.push_back( dset );
  return 0;
}

// Analysis_Hist::setupDimension()
/** Given an ArgList containing name,[min,max,step,bins,col,N], set up a 
  * coordinate with that name and parameters min, max, step, bins.
  * If '*' or not specified, a default value will be set later.
  * \return 1 if error occurs, 0 otherwise.
  */
int Analysis_Hist::setupDimension(ArgList &arglist, DataSet *dset) {
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
      dim.SetMin( DS_Math::Min(*dset) );
    else
      dim.SetMin( default_dim_.Min() );
  }
  // If no max arg and no default max arg, get max from dataset
  if (!maxArg) {
    if (!maxArgSet_)
      dim.SetMax( DS_Math::Max(*dset) );
    else
      dim.SetMax( default_dim_.Max() );
  }
  // Check that min < max
  if (dim.Min() >= dim.Max()) {
    mprinterr("Error: Hist: Dimension %s: min (%lf) must be less than max (%lf).\n",
              dim.c_str(), dim.Min(), dim.Max());
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
  hist_->AddDimension( dim );

  return 0;
}

// Analysis_Hist::Setup()
/** Set up histogram with specified data sets. */
Analysis::RetType Analysis_Hist::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  debug_ = debugIn;

  // Set up histogram DataSet
  hist_ = (Histogram*) datasetlist->AddSet( DataSet::HIST, 
                                            analyzeArgs.GetStringKey("name"), 
                                            "Hist");
  hist_->SetDebug(debug_);
  //hist_ = new Histogram( );

  // Keywords
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

  // Datasets
  // Treat all remaining arguments as dataset names.
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
  // If one or two dimensions and not gnuplot/circular, use DataFileList 
  // for writing.
  if (dimensionArgs_.size() < 3 && !circular_ && !gnuplot_) {
    outfile_ = DFLin->AddDataFile(outfilename_, analyzeArgs);
    if (outfile_ == 0) {
      mprinterr("Error: Could not create output file %s\n", outfilename_.c_str());
      return Analysis::ERR;
    }
    // NOTE: Do not add hist DataSet to outfile here; unlike other DataSets
    // the dimension of hist is not finalized until dimensions are set up in
    // Analyze()
    //outfile_->AddSet(hist_);
  }

  mprintf("\tHist: %s: Set up for %zu dimensions using the following datasets:\n", 
          outfilename_.c_str(), dimensionArgs_.size());
  mprintf("\t      [ ");
  for (std::vector<DataSet*>::iterator ds=histdata_.begin(); ds!=histdata_.end(); ++ds)
    mprintf("%s ",(*ds)->Legend().c_str());
  mprintf("]\n");
  if (calcFreeE_)
    mprintf("\t      Free energy will be calculated from bin populations at %lf K.\n",Temp_);
  if (circular_ || gnuplot_) {
    mprintf("\tWarning: gnuplot and/or circular specified; advanced grace/gnuplot\n");
    mprintf("\t         formatting disabled.\n");
    if (circular_)
      mprintf("\t      circular: Output coordinates will be wrapped.\n");
    if (gnuplot_)
      mprintf("\t      gnuplot: Output will be in gnuplot-readable format.\n");
  }
  if (normalize_)
    mprintf("\t      norm: Bins will be normalized to 1.0.\n");

  return Analysis::OK;
}

// Analysis_Hist::Analyze()
Analysis::RetType Analysis_Hist::Analyze() {
  // Set up dimensions
  // Size of histdata and dimensionArgs should be the same
  for (unsigned int hd = 0; hd < histdata_.size(); hd++) {
    if ( setupDimension(dimensionArgs_[hd], histdata_[hd]) ) 
      return Analysis::ERR;
  }

  // Check that the number of data points in each dimension are equal
  int Ndata = -1;
  for (std::vector<DataSet*>::iterator ds = histdata_.begin(); ds != histdata_.end(); ++ds)
  {
    //mprintf("DEBUG: DS %s size %i\n",histdata[hd]->Name(),histdata[hd]->Xmax()+1);
    if (Ndata==-1)
      Ndata = (*ds)->Size();
    else {
      if (Ndata != (*ds)->Size()) {
        mprinterr("Error: Hist: Dataset %s has inconsistent # data points (%i), expected %i.\n",
                  (*ds)->Legend().c_str(), (*ds)->Size(), Ndata);
        return Analysis::ERR;
      }
    }
  }
  mprintf("\tHist: %i data points in each dimension.\n", Ndata);

  std::vector<double> coord( hist_->NumDimension() );
  for (int n=0; n < Ndata; n++) {
    std::vector<double>::iterator coord_it = coord.begin();
    for (std::vector<DataSet*>::iterator ds = histdata_.begin(); ds != histdata_.end(); ++ds) {
      *coord_it = (*ds)->Dval( n );
      ++coord_it;
    }
    hist_->BinData( coord );
  }

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
