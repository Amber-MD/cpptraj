#include "Histogram.h"
#include <cstdlib>
#include <cstring> //strcpy
#include <cmath> // log
#include "CpptrajStdio.h"
// Histogram

// CONSTRUCTOR
Histogram::Histogram() {
  debug=0;
  Dimension=NULL;
  numDimension=0;
  Bins=NULL;
  currentBin=0;
  isCircular=0;
  BinIndices=NULL;
}

// DESTRUCTOR
Histogram::~Histogram() {
  if (Dimension!=NULL) free(Dimension);
  if (Bins!=NULL) free(Bins);
  if (BinIndices!=NULL) free(BinIndices);
}

// Histogram::SetDebug()
void Histogram::SetDebug(int debugIn) {
  debug = debugIn;
  if (debug>0) mprintf("Histogram DEBUG LEVEL SET TO %i\n", debug);
}

// Histogram::AddDimension()
/** Add a dimension to the histogram with the given min, max, and step.
  * Number of bins will be calculated.
  */
int Histogram::AddDimension(char *labelIn, double minIn, double maxIn, double stepIn) {
  if (debug>0) mprintf("\t\tCalculating bins.\n");
  if (stepIn<=0) {
    mprinterr("Error: Histogram: Dimension %s: step <=0!\n",labelIn);
    return 1;
  }
  double temp = ((maxIn - minIn) / stepIn);
  temp = ceil(temp);
  int dbins = (int) temp;
  return ( AddDimension(labelIn,minIn,maxIn,stepIn,dbins) );
}

// Histogram::AddDimension()
/** Add a dimension to the histogram with the given min, max, step,
  * and number of bins.
  * NOTE: NO ERROR CHECKING IS PERFORMED.
  */
int Histogram::AddDimension(char *labelIn, double minIn, double maxIn, 
                            double stepIn, int binsIn) {
  int offset;

  Dimension = (dimensionType*) realloc( Dimension, (numDimension+1) * sizeof(dimensionType) );
  if (labelIn==NULL)
    Dimension[numDimension].label[0]='\0';
  else
    strcpy(Dimension[numDimension].label, labelIn);
  Dimension[numDimension].min = minIn;
  Dimension[numDimension].max = maxIn;
  Dimension[numDimension].step = stepIn;
  Dimension[numDimension].bins = binsIn;
  if (debug>0) {
    mprintf("\t\t%s: %8.3lf:%8.3lf:%6.2lf:%i\n",Dimension[numDimension].label, 
            Dimension[numDimension].min, Dimension[numDimension].max, 
            Dimension[numDimension].step, Dimension[numDimension].bins);
  }
  numDimension++;
  // NOTE: Should the following be its own routine?
  // Recalculate offsets for all dimensions starting at farthest coord. This 
  // follows column major ordering. 
  offset=1; 
  for (int i = numDimension-1; i >= 0; i--) {
    if (debug>0) mprintf("\tHistogram: %s offset is %i\n",Dimension[i].label, offset);
    Dimension[i].offset=offset;
    offset *= Dimension[i].bins; 
  }
  // offset should now be equal to the total number of bins across all dimensions
  if (debug>0) mprintf("\tHistogram: Total Bins = %i\n",offset);
  numBins = offset;
  // Allocate space for bins
  Bins = (double*) realloc( Bins, numBins * sizeof(double));
  // Set all bins to 0.
  for (int n=0; n < numBins; n++)
    Bins[n]=0; 
  return 0;
}

// Histogram::BinData()
/** Given an array of data with the same dimensionality as the histogram, bin
  * the data.
  */
int Histogram::BinData(double *Data) {
  int index=0,idx;
  double coord;

  // DEBUG
  //mprintf("Binning [");
  //for (int i=0; i<numDimension; i++) mprintf("%lf,",Data[i]);
  //mprintf("]\n");
  // Loop over defined dimensions. 
  // Calculate an index into Bins based on precalcd offsets for dimensions.
  // Populate bin.
  if (debug>1) mprintf("\t{");
  for (int n=0; n<numDimension; n++) {
    // Check if Data is out of bounds for this coordinate
    if (Data[n]>Dimension[n].max || Data[n]<Dimension[n].min) {
      index = -1;
      break;
    }
    // Calculate index (idx) for this particular dimension 
    coord = Data[n] - Dimension[n].min;
    coord = coord / Dimension[n].step;
    //modf(coord,&intpart);
    idx=(int) coord;
    //idx=(int) intpart;
    if (debug>1) mprintf(" [%s:%lf (%i)],",Dimension[n].label,Data[n],idx);

    /* 
    // Check if idx is out of bounds for this dimension 
    if ( (idx<0)||(idx >= Dimension[n].numBins) ) {
      if (debug>1) mprintf("Out of bounds.\n");
      index=-1;
      break;
    }
    */

    // Calculate overall index in Bins, offset has already been calc.
    index+=(idx*Dimension[n].offset);
  }

  // If index was successfully calculated, populate bin 
  if (index>-1 && index < numBins) {
    if (debug>1) mprintf(" |index=%i",index);
    Bins[index]++;
  } else {
    mprintf("\tWarning: Coordinates out of bounds %i {",index);
    for (int n=0; n<numDimension; n++) mprintf("%lf",Data[n]);
    mprintf("}\n");
  }

  if (debug>1) mprintf("}\n");
  return 0;
}

// count2coord()
/** Given an array of integers corresponding to bin indices, calculate the
  * coordinates of those bins and print them out.
  */
void Histogram::count2coord(int *count) {
  for (int i=0; i<numDimension; i++)
    mprintf("%lf ", (count[i]*Dimension[i].step)+Dimension[i].min );
}

// Histogram::PrintBins()
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
  * if SD is not NULL, the standard deviation array will also be printed.
  */
void Histogram::PrintBins(bool circularIn, bool gnuplot) {
  int *count,idx,index,ndim;
  int circular = 0;
  bool loop;
  //int *iBins;
  //double *fBins;
  if (circularIn) circular=1;

  //iBins=NULL; fBins=NULL;

  count=(int*) malloc(numDimension*sizeof(int));
  // Initialize count array according to circular
  if (circular)
    for (int n=0; n<numDimension; n++) count[n]=-1;
  else
    for (int n=0; n<numDimension; n++) count[n]=0;

  if (debug>0) {
    if (circular)
      mprintf("Printing %i bins in circular fashion.\n",numBins);
    else
      mprintf("Printing %i bins.\n",numBins);
  }

  // Set format 
/*  
  switch (binType) {
    case 0 : iBins=(int*) S->Bins;    break;
    case 1 : fBins=(double*) S->Landscape; break;
    default: fprintf(stderr,"binType %i not recognized.\n",binType); return;
  }
*/

  loop=true;
  while (loop) {
    index=0;
    // Calculate index, converting wrapped indices to actual indices 
    for (int n=0; n<numDimension; n++) {
      //mprinterr(" %i",count[n]);
      if (count[n]==-1) idx=Dimension[n].bins-1;
      else if (count[n]==Dimension[n].bins) idx=0;
      else idx=count[n];
      index+=(idx*Dimension[n].offset);
    }
    //mprinterr(" = %i\n",index);

    // If we dont care about zero bins or bin pop > 0, output
    //if (S->nozero==0 || S->Bins[index]!=0) {
      // Output
      count2coord(count);
      //switch (binType) {
        //case 0 :
          mprintf("%lf\n",Bins[index]);
        //break;
        //case 1 :
        //  fprintf(stdout,"%lf",S->Landscape[index]);
        //  if (S->SD!=NULL) fprintf(stdout," %lf",S->SD[index]);
        //  fprintf(stdout,"\n");
        //break;
        //default: fprintf(stderr,"binType %i not recognized.\n",binType); return;
      //}
    //}
/*    
    count2coord(count,S->numCoord,S->Coord);
    if (iBins!=NULL) fprintf(stdout,"%i",iBins[index]);
    if (fBins!=NULL) fprintf(stdout,"%lf",fBins[index]);
    if (S->SD!=NULL) fprintf(stdout," %lf",S->SD[index]);
    fprintf(stdout,"\n");
*/

    // Increment highest order coord
    count[numDimension-1]++;
    // If gnuplot, print extra space when highest order coord cycles
    if (gnuplot) {
      if (count[numDimension-1]==Dimension[numDimension-1].bins+circular) mprintf("\n");
    }
    // Increment other coords if necessary
    for (ndim=numDimension-1; ndim>0; ndim--)
      if (count[ndim]==Dimension[ndim].bins+circular) {count[ndim]=0-circular; count[ndim-1]++;}

    if (count[0]==Dimension[ndim].bins+circular) loop=false;
  }

  free(count);
}

// Histogram::BinStart()
/** Set current bin to 0 and initialize indices. If isCircularIn is true the
  * bin indices will wrap in each dimension.
  */
void Histogram::BinStart(bool isCircularIn) {
  currentBin=0;
  if (isCircularIn) 
    isCircular=1;
  else
    isCircular=0;
  BinIndices = (int*) realloc( BinIndices, numDimension * sizeof(int));
  if (isCircular) 
    for (int dim=0; dim < numDimension; dim++) BinIndices[dim]=-1;
  else
    for (int dim=0; dim < numDimension; dim++) BinIndices[dim]=0;
}

// Histogram::CurrentBinData()
/** Return the data at current bin.  */
double Histogram::CurrentBinData() {
  int index, idx;

  if (isCircular) {
  // If circular, some bins will be accessed more than once so the current 
  // bin must be calculated from indices.
    index=0;
    // Calculate index, converting wrapped indices to actual indices 
    for (int n=0; n<numDimension; n++) {
      //mprintf(" %i",BinIndices[n]);
      if (BinIndices[n]==-1) idx=Dimension[n].bins-1;
      else if (BinIndices[n]==Dimension[n].bins) idx=0;
      else idx=BinIndices[n];
      index+=(idx*Dimension[n].offset);
    }
    //mprintf(" = %i\n",index);
    currentBin = index;
  }
  return Bins[currentBin]; 
}

// Histogram::CurrentBinCoord()
/** Set numDimension coordinates corresponding to the current bin indices.
  */
void Histogram::CurrentBinCoord(double *coord) {
  for (int i=0; i<numDimension; i++)
    coord[i] = (BinIndices[i]*Dimension[i].step)+Dimension[i].min;
}

// Histogram::NextBin()
/** Increment bin and indices. 
  * \return 1 if final bin reached.
  */
int Histogram::NextBin() {
  int dim;

  // Increment bin # - only matters if not circular.
  currentBin++;

  // Increment highest order coord
  BinIndices[numDimension-1]++;
  // Increment other coords if necessary
  for (dim=numDimension-1; dim>0; dim--) {
    if (BinIndices[dim]==Dimension[dim].bins+isCircular) {
      BinIndices[dim]=0-isCircular;
      BinIndices[dim-1]++;
    }
  }

  // If last bin, return 1
  if (BinIndices[0] == Dimension[dim].bins+isCircular) return 1;

  return 0;
}

// Histogram::NBins()
/** Return number of bins for the given dimension */
int Histogram::NBins(int dim) {
  if (dim < 0 || dim >= numDimension) return -1;
  return Dimension[dim].bins;
}

// Histogram::Step() 
/** Return step value for given dimension.  */
double Histogram::Step(int dim) {
  if (dim < 0 || dim >= numDimension) return -1.0;
  return Dimension[dim].step;
}

// Histogram::Min()
/** Return min value for given dimension.  */
double Histogram::Min(int dim) {
if (dim < 0 || dim >= numDimension) return -1.0;
  return Dimension[dim].min;
}

// Histogram::Max() 
/** Return max value for given dimension.  */
double Histogram::Max(int dim) {
if (dim < 0 || dim >= numDimension) return -1.0;
  return Dimension[dim].max;
}

// Histogram::Label()
/** Return label for given dimension.  */
char *Histogram::Label(int dim) {
  if (dim < 0 || dim >= numDimension) return NULL;
  return Dimension[dim].label;
}

// Histogram::BinTotal()
/** Return sum of all bin values */
double Histogram::BinTotal() {
  double sumBin = 0;
  for (int bin=0; bin < numBins; bin++)
    sumBin += Bins[bin];
  mprintf("\tHistogram: Sum of all bins is %lf\n",sumBin);
  return sumBin;
}

// Histogram::CalcFreeE()
/** Calculate free energy based on bin populations.  */
int Histogram::CalcFreeE(double T, int refbin) {
  double KT, binmax, temp, ceiling;

  mprintf("\tHistogram: Calculating free E at %lf K.\n",T);
  KT=(-.001985886596*T);

  // Find most populated bin for G=0
  binmax = Bins[0];
  for (int bin=1; bin < numBins; bin++) 
    if (Bins[bin] > binmax) binmax = Bins[bin];
  mprintf("\t           Bins max is %.0lf\n",binmax);
  if (binmax==0) {
    mprinterr("Histogram: Cannot calc free E, no bins populated!\n");
    return 1;
  }

  // If requested, set up reference bin other than max
  if (refbin>-1) {
    if (Bins[refbin] > 0) {
      binmax=(double)Bins[refbin];
      mprintf("\t           Calculating free E w.r.t bin %i, population %lf\n",refbin,binmax);
    } else
      mprintf("Warning: Reference bin %i has no population. Using %lf\n",refbin,binmax);
  }

  // Set artificial ceiling for bins with 0 population. Make it equivalent
  // to a bin population of 0.5. 
  temp=0.5/binmax;
  ceiling=log(temp)*KT;
  //ceiling*=1.10;
  mprintf("\t           Artificial ceiling is %lf kcal/mol.\n",ceiling);

  // Calculate free E based on populations
  for (int i=0; i<numBins; i++) {
    temp = Bins[i];                // Store Bin population in temp
    if (temp>0) {
      temp/=binmax;                // bin/max
      Bins[i]=log(temp)*KT;   // -R*T*ln(bin/max)
    } else
      Bins[i]=ceiling;        // Artificial ceiling for 0 pop bins.
  }

  return 0;
}

// Histogram::Normalize()
/** Normalize bins so that sum over all bins is 1.0 */
int Histogram::Normalize() {
  double sum = 0.0;
  mprintf("\tHistogram: Normalizing bin populations to 1.0\n");
  for (int i=0; i < numBins; i++)
    sum += Bins[i];
  mprintf("\t           Sum over all bins is %lf\n",sum);
  if (sum == 0.0) {
    mprinterr("Error: Histogram::Normalize: Sum over bin populations is 0.0\n");
    return 1;
  }
  for (int i=0; i < numBins; i++)
    Bins[i] /= sum;
  return 0;
}

