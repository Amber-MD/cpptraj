#include "Analysis_Hist.h"
#include "CpptrajStdio.h"
#include <cstdio> // sprintf
// Analysis_Hist

// CONSTRUCTOR
Hist::Hist() :  
  calcFreeE_(false),
  Temp_(-1.0),
  normalize_(false),
  gnuplot_(false),
  circular_(false),
  outfilename_(NULL),
  minArgSet_(false),
  maxArgSet_(false)
{}

// Hist::CheckDimension()
/** Given an argument with format: DataSet_Name[:min:max:step:bins], check
  * that DataSet_Name exists and is valid. Add the argument to 
  * dimensionArgs and the corresponding dataset to histdata.
  */
int Hist::CheckDimension(char *input, DataSetList *datasetlist) {
  ArgList arglist;
  // Separate input string by ':'
  arglist.SetList(input, ":");
  if (arglist.Nargs()<1) {
    mprintf("Warning: Hist::CheckDimension: No arguments found in input: %s\n",input);
    return 1;
  }

  // First argument should specify dataset name
  if (debug>0) mprintf("\tHist: Setting up histogram dimension using dataset %s\n",
                       arglist.ArgAt(0));
  DataSet *dset = datasetlist->Get( arglist.ArgAt(0) );
  if (dset == NULL) {
    mprintf("\t      Dataset %s not found.\n",arglist.ArgAt(0));
    return 1;
  }

  // Check that dataset is not string
  if (dset->Type()==DataSet::STRING) {
    mprintf("Error: Hist: Cannot histogram dataset %s, type STRING.\n", dset->Name());
    return 1;
  }

  dimensionArgs_.push_back( arglist );
  histdata_.push_back( dset );
  return 0;
}

// Hist::setupDimension()
/** Given an ArgList containing name,[min,max,step,bins,col,N], set up a 
  * coordinate with that name and parameters min, max, step, bins.
  * If '*' or not specified, a default value will be set later.
  * \return 1 if error occurs, 0 otherwise.
  */
int Hist::setupDimension(ArgList &arglist, DataSet *dset) {
  Dimension dim;
  bool minArg = false;
  bool maxArg = false;

  if (debug>1)
    arglist.PrintList();

  // Set up dimension name
  // NOTE: arglist[0] should be same as dset name from CheckDimension 
  dim.SetLabel( arglist[0] );

  // Cycle through coordinate arguments. Any argument left blank will be 
  // assigned a default value later.
  for (int i=1; i<arglist.Nargs(); i++) {
    if (debug>1) mprintf("    DEBUG: setupCoord: Token %i (%s)\n",i,arglist.ArgAt(i));
    // Default explicitly requested
    if (arglist.ArgIs(i,"*")) continue;
    switch (i) {
      case 1 : dim.SetMin( arglist.ArgToDouble(i) ); minArg=true; break;
      case 2 : dim.SetMax( arglist.ArgToDouble(i) ); maxArg=true; break;
      case 3 : dim.SetStep( arglist.ArgToDouble(i) ); break;
      case 4 : dim.SetBins( arglist.ArgToInteger(i) ); break;
    }
  }

  // If no min arg and no default min arg, get min from dataset
  if (!minArg) {
    if (!minArgSet_) 
      dim.SetMin( dset->Min() );
    else
      dim.SetMin( default_dim_.Min() );
  }
  // If no max arg and no default max arg, get max from dataset
  if (!maxArg) {
    if (!maxArgSet_)
      dim.SetMax( dset->Max() );
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
  hist_.AddDimension( dim );

  return 0;
}

// Hist::Setup()
/** Set up histogram with specified data sets.
  * usage: hist(ogram) <dataset_name>[:min:max:step:bins] ...
  *        [free <temperature>] [norm] [gnu] [circular] out <filename>
  *        min <min> max <max> step <step> bins <bins>
  */
int Hist::Setup(DataSetList *datasetlist) {
  char *datasetstring;

  hist_.SetDebug(debug);
  // Keywords
  outfilename_ = analyzeArgs.getKeyString("out",NULL);
  if (outfilename_==NULL) {
    mprintf("Error: Hist: No output filename specified.\n");
    return 1;
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
  while ( (datasetstring = analyzeArgs.getNextString())!=NULL )
    if (CheckDimension( datasetstring,datasetlist )) return 1;
    //if (setupDimension(datasetstring,datasetlist)) return 1;

  mprintf("\tHist: %s: Set up for %zu dimensions using the following datasets:\n", 
          //outfilename, hist.NumDimension());
          outfilename_, dimensionArgs_.size());
  mprintf("\t      [ ");
  for (std::vector<DataSet*>::iterator ds=histdata_.begin(); ds!=histdata_.end(); ++ds)
    mprintf("%s ",(*ds)->Name());
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

  return 0;
}

// Hist::Analyze()
int Hist::Analyze() {
  // Set up dimensions
  // Size of histdata and dimensionArgs should be the same
  for (unsigned int hd = 0; hd < histdata_.size(); hd++) {
    if ( setupDimension(dimensionArgs_[hd], histdata_[hd]) ) 
      return 1;
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
                  (*ds)->Name(), (*ds)->Size(), Ndata);
        return 1;
      }
    }
  }
  mprintf("\tHist: %i data points in each dimension.\n", Ndata);

  std::vector<double> coord( hist_.NumDimension() );
  for (int n=0; n < Ndata; n++) {
    std::vector<double>::iterator coord_it = coord.begin();
    for (std::vector<DataSet*>::iterator ds = histdata_.begin(); ds != histdata_.end(); ++ds) {
      *coord_it = (*ds)->Dval( n );
      ++coord_it;
    }
    hist_.BinData( coord );
  }

  //hist.PrintBins(false,false);
  return 0;
}

// Hist::Print()
/** Convert 1D and 2D histograms to datafiles, otherwise use histogram
  * native output to print.
  */
void Hist::Print(DataFileList *datafilelist) {
  DataFile *outfile=NULL;

  // Calc free energy if requested
  if (calcFreeE_) hist_.CalcFreeE(Temp_, -1);

  // Normalize if requested
  if (normalize_) hist_.Normalize();

  if (hist_.NumDimension() == 1 && !gnuplot_ && !circular_) {
    // For 1 dimension just need to hold bin counts
    hist_.Print_1D( histout_ );
    outfile = datafilelist->Add( outfilename_, histout_.GetDataSetN(0) );
    if (outfile==NULL) {
      mprinterr("Error creating 1D histogram output for file %s\n",outfilename_);
      return;
    }
    outfile->SetXlabel( (char*)hist_[0].c_str() );
    outfile->ProcessArgs("ylabel Count");
    outfile->SetCoordMinStep(hist_[0].Min(), hist_[0].Step(), -1, -1);

  } else if (hist_.NumDimension() == 2 && !gnuplot_ && !circular_) {
    // A DataSet is created for each Y value (dimension 1).
    hist_.Print_2D( histout_ );
    for (DataSetList::const_iterator dset = histout_.begin();
                                     dset != histout_.end(); ++dset)
    {
      outfile = datafilelist->Add( outfilename_, *dset);
    }
    if (outfile==NULL) {
      mprinterr("Error creating %iD histogram output for file %s\n",
                hist_.NumDimension(), outfilename_);
      return;
    }
    outfile->SetXlabel( (char*)hist_[0].c_str() );
    outfile->SetYlabel( (char*)hist_[1].c_str() );
    outfile->SetCoordMinStep(hist_[0].Min(), hist_[0].Step(), hist_[1].Min(), hist_[1].Step());
    outfile->ProcessArgs("noxcol");
    outfile->ProcessArgs("usemap");
    outfile->ProcessArgs("nolabels");
  } else {
    // If > two dimensions, create 1 coord dataset for each dimension plus
    // 1 dataset to hold bin counts.
    //hist_.Print_ND( histout_, circular_ );
    hist_.PrintBins(outfilename_, circular_, gnuplot_);
  }
}
