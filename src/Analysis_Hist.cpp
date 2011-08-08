#include "Analysis_Hist.h"
#include "CpptrajStdio.h"
#include <cmath> // ceil
#include <cstdlib> //atof, atoi
// Analysis_Hist

// CONSTRUCTOR
Hist::Hist() { 
  calcFreeE = false;
  Temp = -1.0;
  normalize = false;
  gnuplot = false;

  min = 0.0;
  max = 0.0;
  step = -1.0;
  bins = -1;

  Ndata = -1;
}

// DESTRUCTOR
Hist::~Hist() {

}

/* Hist::SetupDimension
 * Given a string with format name:min:max:step:bins:col:N, set up a 
 * coordinate with that name and parameters min, max, step, bins.
 * If '*' or not specified, a default value will be set later.
 * Return 1 if input is NULL or error occurs, 0 otherwise.
 */
int Hist::SetupDimension(char *input, DataSetList *datasetlist) {
  ArgList *arglist;
  double dmin,dmax,dstep;
  int dbins;
  DataSet *dset = NULL;

  if (input==NULL) return 1;
 
  // Separate input string by ':'
  arglist = new ArgList(input, ":");
 
  if (arglist->Nargs()<1) {
    mprintf("Warning: Hist::SetupDimension: No arguments found in input: %s\n",input);
    delete arglist;
    return 1;
  }

  // First argument should specify dataset name
  if (debug>0) mprintf("\tHist: Setting up histogram dimension using dataset %s\n",
                       arglist->Arg(0));
  dset = datasetlist->Get(arglist->Arg(0));
  if (dset == NULL) {
    mprintf("\t      Dataset not found.\n");
    delete arglist;
    return 1;
  }

  // Check that the number of data points in each dimension are equal
  if (Ndata==-1)
    Ndata = dset->Xmax();
  else {
    if (Ndata != dset->Xmax()) {
      mprinterr("Error: Hist: Dataset %s does not have correct number of data points (%i).\n",
                dset->Name(),Ndata);
      delete arglist;
      return 1;
    }
  }

  /* 
  // Get range list from first arg 
  range1=getRange(arglist[0],&r1,S->Data,S->Ndata);
  if (range1==NULL) {
    fprintf(stderr,"Error: setupCoord: Could not set up datasets based on argument: %s\n",arglist[0]);
    return;
  }
  */

  // Set up dimension defaults
  //dmin=0.0; dmax=0.0; dstep=-1.0; dbins=-1; 
  dmin = min; dmax = max; dstep = step; dbins = bins;

  // Cycle through coordinate arguments. Any argument left blank will be 
  // assigned a default value later.
  for (int i=1; i<arglist->Nargs(); i++) {
    if (debug>0) mprintf("    DEBUG: setupCoord: Token %i (%s)\n",i,arglist->Arg(i));
    // Default explicitly requested
    if (arglist->ArgIs(i,"*")) continue;
    switch (i) {
      case 1 : dmin = atof(arglist->Arg(i)); break;
      case 2 : dmax = atof(arglist->Arg(i)); break;
      case 3 : dstep= atof(arglist->Arg(i)); break;
      case 4 : dbins= atoi(arglist->Arg(i)); break;
    }
  }
  delete arglist;

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
      if (debug>0) mprintf("    Calculating step.\n");
      if (dbins<=0) {
        mprinterr("Error: Hist: Dimension %s: bins <=0!\n",dset->Name());
        return 1;
      }
      dstep = dmax - dmin;
      dstep = dstep / dbins;
    } else if (dbins==-1) {
      if (debug>0) mprintf("    Calculating bins.\n");
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
 *        [free <temperature>] [norm] [gnu]
 *        min <min> max <max> step <step> bins <bins>
 */
int Hist::Setup(DataSetList *datasetlist) {
  char *datasetstring;

  debug=1;
  // Keywords
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
    SetupDimension(datasetstring,datasetlist);

  mprintf("\tHist: Set up for %i dimensions.\n",hist.NumDimension());

  return 0;
}

/* Hist::Analyze()
 */
int Hist::Analyze() {
  double *coord;

  coord = (double*) malloc( hist.NumDimension() * sizeof(double));

  for (int n=0; n < Ndata; n++) {
    for (int hd=0; hd < (int)histdata.size(); hd++)
      coord[hd] = histdata[hd]->Get((void*)(coord+hd), n);
    hist.BinData(coord);
  }

  free(coord);

  hist.PrintBins(0,false);
  return 0;
}
