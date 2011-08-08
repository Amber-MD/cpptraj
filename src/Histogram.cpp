#include "Histogram.h"
#include <cstdlib>
#include <cstring> //strcpy
#include "CpptrajStdio.h"
// Histogram

// CONSTRUCTOR
Histogram::Histogram() {
  debug=0;
  Dimension=NULL;
  numDimension=0;
  Bins=NULL;
}

// DESTRUCTOR
Histogram::~Histogram() {
  if (Dimension!=NULL) free(Dimension);
  if (Bins!=NULL) free(Bins);
}

/* Histogram::SetDebug()
 */
void Histogram::SetDebug(int debugIn) {
  debug = debugIn;
  if (debug>0) mprintf("HISTOGRAM DEBUG LEVEL SET TO %i\n", debug);
}

/* Histogram::AddDimension()
 * Add a dimension to the histogram with the given min, max, step,
 * and number of bins.
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
  mprintf("\t\t%s: %8.3lf:%8.3lf:%6.2lf:%i\n",Dimension[numDimension].label, 
          Dimension[numDimension].min, Dimension[numDimension].max, 
          Dimension[numDimension].step, Dimension[numDimension].bins);
  numDimension++;
  // NOTE: Should the following be its own routine?
  // Recalculate offsets for all dimensions starting at farthest coord. This 
  // follows column major ordering. 
  offset=1; 
  for (int i = numDimension-1; i >= 0; i--) {
    if (debug>0) mprintf("\tHIST: %s offset is %i\n",Dimension[i].label, offset);
    Dimension[i].offset=offset;
    offset *= Dimension[i].bins; 
  }
  // offset should now be equal to the total number of bins across all dimensions
  mprintf("\tHIST: Total Bins = %i\n",offset);
  numBins = offset;
  // Allocate space for bins
  Bins = (int*) realloc( Bins, numBins * sizeof(int));
  // Set all bins to 0.
  for (int n=0; n < numBins; n++)
    Bins[n]=0; 
  return 0;
}

/* Histogram::BinData()
 * Given an array of data with the same dimensionality as the histogram, bin
 * the data.
 */
int Histogram::BinData(double *Data) {
  int index=0,idx;
  double coord;

  // Loop over defined dimensions. 
  // Calculate an index into Bins based on precalcd offsets for dimensions.
  // Populate bin.
  if (debug>1) mprintf("{");
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
    // Check if i is out of bounds for this dimension 
    if ( (i<0)||(i >= Dimension[n].numBins) ) {
      if (debug>1) fprintf(stdout,"Out of bounds.\n");
      index=-1;
      break;
    }
    */

    // Calculate overall index in Bins, offset has already been calc.
    index+=(idx*Dimension[n].offset);
  }

  // If index was successfully calculated, populate bin */
  if (index!=-1) {
    if (debug>1) mprintf(" |index=%i",index);
    Bins[index]++;
  }

  if (debug>1) mprintf("}\n");
  return 0;
}

/* count2coord()
 * Given an array of integers corresponding to bin indices, calculate the
 * coordinates of those bins and print them out.
 */
void Histogram::count2coord(int *count) {
  for (int i=0; i<numDimension; i++)
    mprintf("%lf ", (count[i]*Dimension[i].step)+Dimension[i].min );
}

/* Histogram::PrintBins()
 * This routine will print out the contents of the Bin array in column-major 
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
void Histogram::PrintBins(int circular) {
  int *count,idx,index,ndim;
  bool loop;
  //int *iBins;
  //double *fBins;

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
      //fprintf(stdout,"%i ",count[n]);
      if (count[n]==-1) idx=Dimension[n].bins-1;
      else if (count[n]==Dimension[n].bins) idx=0;
      else idx=count[n];
      index+=(idx*Dimension[n].offset);
    }
    //fprintf(stdout," = %i\n",index);

    // If we dont care about zero bins or bin pop > 0, output
    //if (S->nozero==0 || S->Bins[index]!=0) {
      // Output
      count2coord(count);
      //switch (binType) {
        //case 0 :
          mprintf("%i\n",Bins[index]);
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
    //if (S->gnuplot==1) {
    //  if (count[numDimension-1]==Dimension[numDimension-1].bins+circular) fprintf(stdout,"\n");
    //}
    // Increment other coords if necessary
    for (ndim=numDimension-1; ndim>0; ndim--)
      if (count[ndim]==Dimension[ndim].bins+circular) {count[ndim]=0-circular; count[ndim-1]++;}

    if (count[0]==Dimension[ndim].bins+circular) loop=false;
  }

  free(count);
}

