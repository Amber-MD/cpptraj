#include <cmath>
#include "Analysis_Corr.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Corr::Corr() {
  D1=NULL;
  D2=NULL;
  lagmax=0;
  Nelements=0;
  outfilename=NULL;
}

// Corr::Setup()
/** Expected call: corr <outfilename> <Dataset1> <Dataset2> [lagmax <lagmax>]
  */
int Corr::Setup(DataSetList *datasetlist) {
  // Keywords
  lagmax = analyzeArgs.getKeyInt("lagmax",-1);
  outfilename = analyzeArgs.getKeyString("out",NULL);
  if (outfilename==NULL) {
    mprinterr("Error: Corr: No output filename specified ('out' <filename>).\n");
    return 1;
  }
  
  // Datasets
  char *D1name = analyzeArgs.getNextString();
  if (D1name==NULL) {
    mprinterr("Error: Corr: Must specify 2 dataset names.\n");
    return 1;
  }
  char *D2name = analyzeArgs.getNextString();
  if (D2name==NULL) {
    mprinterr("Error: Corr: Must specify 2 dataset names.\n");
    return 1;
  }
  D1 = datasetlist->Get(D1name);
  if (D1==NULL) {
    mprinterr("Error: Corr: Could not get dataset named %s\n",D1name);
    return 1;
  }
  D2 = datasetlist->Get(D1name);
  if (D2==NULL) {
    mprinterr("Error: Corr: Could not get dataset named %s\n",D2name);
    return 1;
  }

  // Check that D1 and D2 have same # data points.
  // NOTE: Should also check type, this should be performed later!
  Nelements = D1->Capacity(); 
  if (Nelements != D2->Capacity()) {
    mprinterr("Error: Corr: # elements in dataset %s (%i) not equal to\n",D1name,Nelements);
    mprinterr("             # elements in dataset %i (%i)\n",D2name, D2->Capacity());
    return 1;
  }
  if (lagmax==-1) lagmax = Nelements;

  // Setup output dataset
  if (Ct.Setup((char*)"Corr", lagmax)) return 1;

  mprintf("    CORR: Correlation between set %s and set %s, %i elements, max lag %i\n",
          D1name, D2name, Nelements, lagmax);
  mprintf("          Output to %s\n",outfilename);

  return 0;
}

// Corr::Analyze()
int Corr::Analyze() {
  double d1;
  double d2;
  double ct;
  // Calculate averages
  double avg1 = D1->Avg(NULL);
  double avg2 = D2->Avg(NULL);
  //mprintf("Avg1=%lf  Avg2=%lf\n",avg1,avg2);
  
  // Compute normalization
  double sumdiff1_2 = 0;
  double sumdiff2_2 = 0;
  for (int i = 0; i < Nelements; i++) {
    d1 = D1->Dval(i);
    d2 = D2->Dval(i);
    double diff1 = d1 - avg1;
    double diff2 = d2 - avg2;
    sumdiff1_2 += (diff1 * diff1);
    sumdiff2_2 += (diff2 * diff2);
  }
  double norm = sumdiff1_2 * sumdiff2_2;
  if (norm <= 0) {
    mprinterr("Error: Corr: Normalization sqrt <= 0.\n");
    return 1;
  }
  norm = sqrt( norm );

  // Calculate correlation
  for (int lag = 0; lag < lagmax; lag++) {
    ct = 0; 
    for (int j = 0; j < Nelements - lag; j++) {
      d1 = D1->Dval(j);
      d2 = D2->Dval(j+lag);
      ct += ((d1 - avg1) * (d2 - avg2));
    }
    ct /= norm;
    Ct.Add(lag, &ct);
  }

  return 0;
}

// Corr::Print()
void Corr::Print(DataFileList *dfl) {
  dfl->Add(outfilename, &Ct);
}

