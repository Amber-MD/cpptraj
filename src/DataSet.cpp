// DataSet
#include <cstdlib>
#include <cstring>
#include "DataSet.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataSet::DataSet() {
  //fprintf(stderr,"DataSet Constructor.\n");
  name=NULL;
  N=0;
  isDynamic=0;
  current=0;
}

// DESTRUCTOR
DataSet::~DataSet() {
  //fprintf(stderr,"DataSet Destructor\n");
  if (name!=NULL) free(name);
}

/* 
 * DataSet::Setup()
 * Set up common to all data sets. The dataset name should be unique and is
 * checked for in DataSetList prior to this call. Nin is the expected size 
 * of the dataset. If Nin<=0 the dataset will be allocated dynamically.
 */
int DataSet::Setup(char *nameIn, int Nin) {
  // Dataset name
  if (nameIn==NULL) {
    mprintf("Dataset has no name.\n");
    return 1;
  }
  name=(char*) malloc( (strlen(nameIn)+1) * sizeof(char));
  strcpy(name, nameIn);
  // Dataset memory
  N=Nin;
  if (N<=0) {
    isDynamic=1;
    N=0;
  }
  if ( this->Allocate() ) return 1;
  return 0;
}

/*
 * DataSet::Info()
 * Print dataset information.
 */
void DataSet::Info() {
  mprintf("    Data set %s",name);
  mprintf(", size is ");
  if (isDynamic)
    mprintf("dynamic");
  else
    mprintf("%i",N);
  mprintf(", current is %i",current);
  mprintf(".\n");
}

/*
 * DataSet::CheckSet()
 * Return 1 if current==0, which indicates set has not been written to.
 * Otherwise return 0;
 */
int DataSet::CheckSet() {
  if (current==0) return 1;
  return 0;
}

