#include "Analysis_Hist.h"
#include "CpptrajStdio.h"
#include <cmath> // ceil
#include <cstdio> // sprintf
// Analysis_Hist

// CONSTRUCTOR
Hist::Hist() { 
  calcFreeE = false;
  Temp = -1.0;
  normalize = false;
  gnuplot = false;
  circular = false;
  Ndata = -1;
  outfilename=NULL;

  defaultMinSet = false;
  defaultMaxSet = false;
  
  min = 0.0;
  max = 0.0;
  step = -1.0;
  bins = -1;
}

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
  DataSet *dset = datasetlist->Get(arglist.ArgAt(0));
  if (dset == NULL) {
    mprintf("\t      Dataset %s not found.\n",arglist.ArgAt(0));
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
    return 1;
  }

  dimensionArgs.AddArg( input );
  histdata.push_back(dset);
  return 0;
}

// Hist::setupDimension()
/** Given a string with format name:min:max:step:bins:col:N, set up a 
  * coordinate with that name and parameters min, max, step, bins.
  * If '*' or not specified, a default value will be set later.
  * \return 1 if error occurs, 0 otherwise.
  */
int Hist::setupDimension(char *input, DataSet *dset) {
  ArgList arglist;
  double dmin,dmax,dstep;
  bool minArg=false;
  bool maxArg=false;
  int dbins;

  // Separate input string by ':'
  arglist.SetList(input, ":");
  if (arglist.Nargs()<1) {
    mprintf("Warning: Hist::setupDimension: No arguments found in input: %s\n",input);
    return 1;
  }

  // Set up dimension defaults. If any arguments are specified then values
  // we be recalculated.
  if (arglist.Nargs() > 1) {
    dmin=0.0; dmax=0.0; dstep=-1.0; dbins=-1;
  } else {
    dmin = min; dmax = max; dstep = step; dbins = bins;
  }

  // Cycle through coordinate arguments. Any argument left blank will be 
  // assigned a default value later.
  for (int i=1; i<arglist.Nargs(); i++) {
    if (debug>1) mprintf("    DEBUG: setupCoord: Token %i (%s)\n",i,arglist.ArgAt(i));
    // Default explicitly requested
    if (arglist.ArgIs(i,"*")) continue;
    switch (i) {
      case 1 : dmin = arglist.ArgToDouble(i); minArg=true; break;
      case 2 : dmax = arglist.ArgToDouble(i); maxArg=true; break;
      case 3 : dstep= arglist.ArgToDouble(i); break;
      case 4 : dbins= arglist.ArgToInteger(i); break;
    }
  }

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

    // If no min arg and no default min arg, get min from dataset
    if (!minArg && dmin == 0) 
      dmin = dset->Min();
    // If no max arg and no default max arg, get max from dataset
    if (!maxArg && dmax == 0)
      dmax = dset->Max();

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
    //histdata.push_back(dset);

  //}

  //safe_free(range1);
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

  hist.SetDebug(debug);
  // Keywords
  outfilename = analyzeArgs.getKeyString("out",NULL);
  if (outfilename==NULL) {
    mprintf("Error: Hist: No output filename specified.\n");
    return 1;
  }
  Temp = analyzeArgs.getKeyDouble("free",-1.0);
  if (Temp!=-1.0) calcFreeE = true;
  if (analyzeArgs.hasKey("gnu")) gnuplot = true;
  if (analyzeArgs.hasKey("norm")) normalize = true;
  //if (analyzeArgs.hasKey("circular")) circular = true;
  // NOTE: The following may only need to be local
  if (analyzeArgs.Contains("min")) {
    min = analyzeArgs.getKeyDouble("min",0.0);
    defaultMinSet = true;
  }
  if (analyzeArgs.Contains("max")) {
    max = analyzeArgs.getKeyDouble("max",0.0);
    defaultMaxSet = true;
  }
  step = analyzeArgs.getKeyDouble("step",-1.0);
  bins = analyzeArgs.getKeyInt("bins",-1);

  // Datasets
  // Treat all remaining arguments as dataset names. 
  while ( (datasetstring = analyzeArgs.getNextString())!=NULL )
    if (CheckDimension( datasetstring,datasetlist )) return 1;
    //if (setupDimension(datasetstring,datasetlist)) return 1;

  mprintf("\tHist: %s: Set up for %i dimensions using the following datasets:\n", 
          //outfilename, hist.NumDimension());
          outfilename, dimensionArgs.Nargs());
  dimensionArgs.PrintList();
  //mprintf("\t      [ ");
  //for (std::vector<DataSet*>::iterator ds=histdata.begin(); ds!=histdata.end(); ds++)
  //  mprintf("%s ",(*ds)->Name());
  //mprintf("]\n");
  if (calcFreeE)
    mprintf("\t      Free energy will be calculated from bin populations at %lf K.\n",Temp);
  if (circular)
    mprintf("\t      circular: Output coordinates will be wrapped.\n");
  if (normalize)
    mprintf("\t      norm: Bins will be normalized to 1.0.\n");

  return 0;
}

// Hist::Analyze()
int Hist::Analyze() {
  double *coord;
  char *datasetstring;
 
  // Set up dimensions
  for (int hd = 0; hd < (int)histdata.size(); hd++) {
    datasetstring = dimensionArgs.ArgAt(hd);
    if (setupDimension(datasetstring,histdata[hd])) return 1;
  }

  // Check that the number of data points in each dimension are equal
  for (int hd=0; hd < (int)histdata.size(); hd++) {
    //mprintf("DEBUG: DS %s size %i\n",histdata[hd]->Name(),histdata[hd]->Xmax()+1);
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

  coord = new double[ hist.NumDimension() ];
  for (int n=0; n < Ndata; n++) {
    for (int hd=0; hd < (int)histdata.size(); hd++) {
      histdata[hd]->Get((void*)(coord+hd), n);
    }
    hist.BinData(coord);
  }
  delete[] coord;

  //hist.PrintBins(false,false);
  return 0;
}

// Hist::Print()
/** Convert 1D and 2D histograms to datafiles, otherwise use histogram
  * native output to print.
  */
void Hist::Print(DataFileList *datafilelist) {
  DataFile *outfile=NULL;
  double *coord;
  int dim, bin;
  double N;
  bool histloop = true;

  //hist.Info();

  // Calc free energy if requested
  if (calcFreeE) hist.CalcFreeE(Temp,-1);

  // Normalize if requested
  if (normalize) hist.Normalize();

  coord = new double[ hist.NumDimension() ];
  hist.BinStart(circular);

  // For 1 dimension just need to hold bin counts
  if (hist.NumDimension() == 1) {
    outfile = datafilelist->Add(outfilename, histout.Add( DOUBLE, hist.Label(0), "Hist" ));
    bin = 0;
    while (histloop) {
      N = hist.CurrentBinData();
      histout.AddData( bin, &N, 0 );
      bin++;
      if (hist.NextBin()) histloop=false;
    }
    outfile->SetXlabel(hist.Label(0));
    outfile->SetYlabel((char*)"Count");
    outfile->SetCoordMinStep(hist.Min(0),hist.Step(0),hist.Min(1),hist.Step(1));

  // The way that datafile understands 2D data currently:
  //   frame0 set0
  //   frame0 set1
  //   frame0 set2
  //   frame1 set0
  //   ...
  // So need a set for each Y value (dimension 1).
  } else if (hist.NumDimension() == 2) {
    char temp[32];
    for (bin = 0; bin < hist.NBins(1); bin++) {
      sprintf(temp,"%8.3lf",(bin*hist.Step(1))+hist.Min(1));
      outfile = datafilelist->Add(outfilename, histout.AddMultiN( DOUBLE, "", temp, bin ));
    }
    bin = 0; // y coord index
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
      //mprintf("%lf %lf %lf [dim=%i bin=%i]\n",coord[0],coord[1],N,dim,bin);
      histout.AddData( dim, &N, bin);
      bin++;
      if (hist.NextBin()) histloop=false;
    }
    outfile->SetXlabel(hist.Label(0));
    outfile->SetYlabel(hist.Label(1));
    outfile->SetCoordMinStep(hist.Min(0),hist.Step(0),hist.Min(1),hist.Step(1));

  // If > two dimensions, create 1 coord dataset for each dimension plus
  // 1 dataset to hold bin counts.
  } else {
   for (dim = 0; dim < hist.NumDimension(); dim++) {
      outfile = datafilelist->Add(outfilename, 
                                  histout.Add( DOUBLE, histdata[dim]->Name(), "Hist" ));
    }
    outfile = datafilelist->Add(outfilename, histout.Add( DOUBLE, NULL, "Count"));
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
  }

  if (hist.NumDimension() > 1) {
    outfile->SetNoXcol();
    outfile->SetMap();
    outfile->SetNoLabels();
  }
  //hist.PrintBins(circular,false);
  delete []coord;
}
