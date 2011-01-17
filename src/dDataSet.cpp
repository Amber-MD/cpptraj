// dDataSet
#include <cstdio> // sprintf
#include <cstdlib>
#include "dDataSet.h"
#include "PtrajMpi.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
dDataSet::dDataSet() {
//  fprintf(stderr,"DDATASET CONSTRUCTOR\n");
  Data=NULL;
  Xvalues=NULL;
}

// DESTRUCTOR
dDataSet::~dDataSet() {
  //fprintf(stderr,"dDataSet Destructor\n");
  if (Data!=NULL) free(Data);
  if (Xvalues!=NULL) free(Xvalues);
}

/* 
 * dDataSet::Allocate()
 * Allocate dataset memory. isDynamic and N are set by DataSet::InternalSetup 
 */
int dDataSet::Allocate() {

  if (isDynamic==0) {
    Data = (double*) calloc( N , sizeof(double));
    Xvalues = (int*) calloc( N , sizeof(int));
  }
  
  return 0;
}

/*
 * dDataSet::Add()
 * Add value to DataSet
 */
void dDataSet::Add( int frame, void *vIn ) {
  double *value;

  value=(double*) vIn;
  if (isDynamic) {
    N++;
    Data=(double*) realloc(Data, N * sizeof(double));
    Xvalues=(int*) realloc(Xvalues, N * sizeof(int));
  } else
    current=frame;

  if (current>=N) {
    mprintf("Error: DataSet %s attempting to write %i, %i allocated.\n",
            name,current,N);
    return;
  }

  Data[current]=*value;
  //Data[frame]=*value;
  Xvalues[current]=frame;
  current++;
  // Increment current. Only matters if isDynamic.
}

/* dDataSet::Write()
 * Formatted write out data for frame to buffer.
 */
void dDataSet::Write(char *buffer, int frame) {
  //if (frame>=N) fprintf(outfile," %12s","NoData"); 
  //fprintf(outfile," %12.4lf",Data[frame]);
  if (frame>=N) 
    sprintf(buffer," %12s","NoData");
  else 
    sprintf(buffer," %12.4lf",Data[frame]);
}

/* dDataSet::Sync()
 * MPI only - sum up this data set across all threads.
 */
int dDataSet::Sync() {
  double *recvBuffer;

  //return 0; // DEBUG
  
  if (worldsize==1) return 0;
#ifdef DEBUG
  dbgprintf("\tSyncing dataset %s\n",name);
#endif

  recvBuffer=(double*) calloc(N, sizeof(double));

  if ( parallel_sum(recvBuffer, Data, N) ) {
    free(recvBuffer);
    return 1;
  }

  // On master, replace Data with recvBuffer
  // NOTE: Only allocate recvBuffer on master?
  if (worldrank==0) {
    free(Data); 
    Data=recvBuffer;
  } else
    free(recvBuffer);

  return 0;
}
