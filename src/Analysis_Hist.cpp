#include "Analysis_Hist.h"
#include "CpptrajStdio.h"
#include <cmath> // ceil
#include <cstdlib> //atof, atoi
#include <cstdio> // sprintf
// Analysis_Hist

// CONSTRUCTOR
Hist::Hist() { 
  calcFreeE = false;
  Temp = -1.0;
  normalize = false;
  gnuplot = false;
  Ndata = -1;
  outfilename=NULL;
  
  min = 0.0;
  max = 0.0;
  step = -1.0;
  bins = -1;
}

// DESTRUCTOR
Hist::~Hist() {

}

/* Hist::setupDimension()
 * Given a string with format name:min:max:step:bins:col:N, set up a 
 * coordinate with that name and parameters min, max, step, bins.
 * If '*' or not specified, a default value will be set later.
 * Return 1 if error occurs, 0 otherwise.
 */
int Hist::setupDimension(char *input, DataSetList *datasetlist) {
  ArgList arglist(input, ":");
  double dmin,dmax,dstep;
  int dbins;
  DataSet *dset = NULL;

  // Separate input string by ':'
  //arglist = new ArgList(input, ":");
 
  if (arglist.Nargs()<1) {
    mprintf("Warning: Hist::setupDimension: No arguments found in input: %s\n",input);
    //delete arglist;
    return 1;
  }

  // First argument should specify dataset name
  if (debug>0) mprintf("\tHist: Setting up histogram dimension using dataset %s\n",
                       arglist.Arg(0));
  dset = datasetlist->Get(arglist.Arg(0));
  if (dset == NULL) {
    mprintf("\t      Dataset not found.\n");
    //delete arglist;
    return 1;
  }

  /* 
  // Get range list from first arg 
  range1=getRange(arglist[0],&r1,S->Data,S->Ndata);
  if (range1==NULL) {
    fprintf(stderr,"Error: setupCoord: Could not set up datasets based on argument: %s\n",arglist[0]);
    return;
  }
  */

  // Check that dataset is not string
  if (dset->Type()==STRING) {
    mprintf("Error: Hist: Cannot histogram dataset %s, type STRING.\n", dset->Name());
    //delete arglist;
    return 1;
  }

  // Set up dimension defaults
  //dmin=0.0; dmax=0.0; dstep=-1.0; dbins=-1; 
  dmin = min; dmax = max; dstep = step; dbins = bins;

  // Cycle through coordinate arguments. Any argument left blank will be 
  // assigned a default value later.
  for (int i=1; i<arglist.Nargs(); i++) {
    if (debug>0) mprintf("    DEBUG: setupCoord: Token %i (%s)\n",i,arglist.Arg(i));
    // Default explicitly requested
    if (arglist.ArgIs(i,"*")) continue;
    switch (i) {
      case 1 : dmin = atof(arglist.Arg(i)); break;
      case 2 : dmax = atof(arglist.Arg(i)); break;
      case 3 : dstep= atof(arglist.Arg(i)); break;
      case 4 : dbins= atoi(arglist.Arg(i)); break;
    }
  }
  //delete arglist;

  // For each dataset index in range1, set up a dimension 
  //for (i=0; i<r1; i++) {
    //if (debug>0) {
    //  fprintf(stderr,"    Coord %i: ",n);
    //  if (min!=NULL) fprintf(stderr,"%lf->",*min); else fprintf(stderr,"*->");
    //  if (max!=NULL) fprintf(stderr,"%lf,",*max); else fprintf(stderr,"*,");
    //  if (step!=-1) fprintf(stderr," step %lf,",step); else fprintf(stderr," step *,");
    //  if (bins!=-1) fprintf(stderr," %i bins.",bins); else fprintf(stderr," * bins.");
    //  fprintf(stderr," Using dataset %s\n",S->Data[column].label);
    //}

    // Check that min < max
    if (dmin >= dmax) {
      mprinterr("Error: Hist: Dimension %s: min (%lf) must be less than max (%lf).\n",dset->Name(),
                dmin,dmax);
      return 1;
    }

    // Attempt to set up bins or step. If both have been defined, use the
    // specified bins and calculate a new step. If neither have been defined,
    // use default bins and calculate step.
    // When calculating bins from a stepsize, round up.
    if (dbins!=-1 && dstep!=-1) {
      mprintf("\tHist: Bins and step have been specified. Recalculating step.\n");
      dstep=-1;
    }

    if (dbins==-1 && dstep==-1) {
      mprintf("\tHist: Bins and step undefined.\n");
      return 1;
    }

    if (dstep==-1) {
      if (debug>0) mprintf("\t\tCalculating step.\n");
      if (dbins<=0) {
        mprinterr("Error: Hist: Dimension %s: bins <=0!\n",dset->Name());
        return 1;
      }
      dstep = dmax - dmin;
      dstep = dstep / dbins;
    } else if (dbins==-1) {
      if (debug>0) mprintf("\t\tCalculating bins.\n");
      if (dstep<=0) {
        mprinterr("Error: Hist: Dimension %s: step <=0!\n",dset->Name());
        return 1;
      }
      double temp = ((dmax - dmin) / dstep);
      temp = ceil(temp);
      dbins = (int) temp;
    }

    mprintf("\tHist: Dim %s: %lf->%lf, step %lf, %i bins.\n", dset->Name(),
            dmin,dmax,dstep,dbins);
    hist.AddDimension(dset->Name(),dmin,dmax,dstep,dbins);
    histdata.push_back(dset);

  //}

  //safe_free(range1);
  return 0;
}

/* Hist::Setup()
 * Set up histogram with specified data sets.
 * usage: hist(ogram) <dataset_name>[:min:max:step:bins] ...
 *        [free <temperature>] [norm] [gnu] out <filename>
 *        min <min> max <max> step <step> bins <bins>
 */
int Hist::Setup(DataSetList *datasetlist) {
  char *datasetstring;

  debug=1;
  // Keywords
  outfilename = analyzeArg->getKeyString("out",NULL);
  if (outfilename==NULL) {
    mprintf("Error: Hist: No output filename specified.\n");
    return 1;
  }
  Temp = analyzeArg->getKeyDouble("free",-1.0);
  if (Temp!=-1.0) calcFreeE = true;
  if (analyzeArg->hasKey("gnu")) gnuplot = true;
  if (analyzeArg->hasKey("norm")) normalize = true;
  // NOTE: The following may only need to be local
  min = analyzeArg->getKeyDouble("min",0.0);
  max = analyzeArg->getKeyDouble("max",0.0);
  step = analyzeArg->getKeyDouble("step",-1.0);
  bins = analyzeArg->getKeyInt("bins",-1);

  // Datasets
  // Treat all remaining arguments as dataset names. 
  while ( (datasetstring = analyzeArg->getNextString())!=NULL )
    if (setupDimension(datasetstring,datasetlist)) return 1;

  mprintf("\tHist: %s: Set up for %i dimensions using the following datasets:\n", 
          outfilename, hist.NumDimension());
  mprintf("\t      [ ");
  for (std::vector<DataSet*>::iterator ds=histdata.begin(); ds!=histdata.end(); ds++)
    mprintf("%s ",(*ds)->Name());
  mprintf("]\n");

  return 0;
}

/* Hist::Analyze()
 */
int Hist::Analyze() {
  double *coord;

  // Check that the number of data points in each dimension are equal
  for (int hd=0; hd < (int)histdata.size(); hd++) {
    mprintf("DEBUG: DS %s size %i\n",histdata[hd]->Name(),histdata[hd]->Xmax()+1);
    if (Ndata==-1)
      Ndata = histdata[hd]->Xmax()+1;
    else {
      if (Ndata != histdata[hd]->Xmax()+1) {
        mprinterr("Error: Hist: Dataset %s has inconsistent # data points (%i), expected %i.\n",
                  histdata[hd]->Name(),histdata[hd]->Xmax()+1,Ndata);
        return 1;
      }
    }
  }
  mprintf("\tHist: %i data points in each dimension.\n",Ndata);

  hist.SetDebug(2);
  coord = (double*) malloc( hist.NumDimension() * sizeof(double));
  for (int n=0; n < Ndata; n++) {
    for (int hd=0; hd < (int)histdata.size(); hd++) {
      histdata[hd]->Get((void*)(coord+hd), n);
    }
    hist.BinData(coord);
  }
  free(coord);

  //hist.PrintBins(false,false);
  return 0;
}

/* Hist::Print()
 * Convert 1D and 2D histograms to datafiles, otherwise use histogram
 * native output to print.
 */
void Hist::Print(DataFileList *datafilelist) {
  DataFile *outfile;
  double *coord;
  int N, dim, bin;
  bool histloop = true;

  hist.Info();

  coord = (double*) malloc(hist.NumDimension() * sizeof(double));
  hist.BinStart(false);

  // If not two dimensions, create 1 coord dataset for each dimension plus
  // 1 dataset to hold bin counts.
  if (hist.NumDimension() != 2) {
    for (dim = 0; dim < hist.NumDimension(); dim++) {
      outfile = datafilelist->Add(outfilename, 
                                  histout.Add( DOUBLE, histdata[dim]->Name(), "Hist" ));
    }
    outfile = datafilelist->Add(outfilename, histout.Add( INT, NULL, "Count"));
    bin = 0;
    while (histloop) {
      hist.CurrentBinCoord(coord);
      N = hist.CurrentBinData();
      for (dim=0; dim<hist.NumDimension(); dim++)
        histout.AddData( bin, coord + dim, dim ); 
        //mprintf("%12.4lf ",coord[dim]);
      histout.AddData( bin, &N, dim );
      bin++;
      //mprintf("%8i\n",N);
      if (hist.NextBin()) histloop=false;
    }

  // The way that datafile understands 2D data currently, each X coord block:
  //   X0:Y0 X0:Y1 X0:Y2 ... X1:Y0 X1:Y1 X1:Y2 ...
  // is stored in its own data set.
  } else {
    char temp[32];
    for (dim = 0; dim < hist.NumBins1D(); dim++) {
      sprintf(temp,"%i",dim);
      outfile = datafilelist->Add(outfilename, histout.AddIdx( INT, temp, dim ));
    }
    bin = 0;
    dim = 0; // x coord index
    hist.CurrentBinCoord(coord);
    double highestcoord = coord[0];
    while (histloop) {
      hist.CurrentBinCoord(coord);
      if (coord[0]!=highestcoord) {
        dim++;
        highestcoord = coord[0];
        bin = 0;
      }
      N = hist.CurrentBinData();
      histout.AddData( bin, &N, dim);
      bin++;
      if (hist.NextBin()) histloop=false;
    }
  }

  outfile->SetNoXcol();
  outfile->SetCoordMinStep(hist.Min(0),hist.Step(0),hist.Min(1),hist.Step(1));
  outfile->SetXlabel(hist.Label(0));
  outfile->SetYlabel(hist.Label(1));
  outfile->SetMap();
  outfile->SetNoLabels();
  //hist.PrintBins(false,false);
  free(coord);
}
