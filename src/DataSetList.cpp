// DataSetList
#include <cstdio> // sprintf
#include <cstdlib>
#include <cstring>
// This also includes basic DataSet class and dataType
#include "DataSetList.h"
#include "CpptrajStdio.h"
// Data types go here
//#include "dDataSet.h"
#include "mapDataSet.h"
#include "stringDataSet.h"
#include "intDataSet.h"

// CONSTRUCTOR
DataSetList::DataSetList() {
  //fprintf(stderr,"DSL Constructor\n");
  DataList=NULL;
  Ndata=0;
  maxFrames=0;
}

// DESTRUCTOR
DataSetList::~DataSetList() {
  int i;
  //fprintf(stderr,"DSL Destructor\n"); 
  if (DataList!=NULL) {
    for (i=0; i<Ndata; i++)
      delete DataList[i];
    free(DataList);
  }
}

/*
 * DataSetList::SetMax()
 * Set the max number frames expected to be read in. Used to preallocate
 * data set sizes in the list.
 */
void DataSetList::SetMax(int expectedMax) {
  maxFrames = expectedMax;
  if (maxFrames<0) maxFrames=0;
}

/* 
 * DataSetList::Get()
 * Return pointer to DataSet with given name
 */
DataSet *DataSetList::Get(char *nameIn) {
  int i;
  for (i=0; i<Ndata; i++)
    if ( strcmp(DataList[i]->Name(),nameIn)==0 ) return DataList[i];
  return NULL;
}

/* 
 * DataSetList::Add()
 * Add a DataSet of specified type, set it up and return pointer to it. 
 * Name the dataset nameIn if specified, otherwise give it a default
 * name based on the given defaultName and dataset #. This routine
 * MUST be called with a default name.
 */ 
DataSet *DataSetList::Add(dataType inType, char *nameIn, const char *defaultName) {
  DataSet *D;
  char tempName[32];

  D=NULL;

  // Require all calls provide a default name
  if (defaultName==NULL) {
    mprinterr("Internal Error: DataSetList::Add() called without default name.\n");
    return NULL;
  }

  // If no name given make one up based on the given default 
  if (nameIn==NULL) {
    sprintf(tempName,"%s_%05i",defaultName,Ndata);
    nameIn=tempName;
  }
  // Check if dataset name is already in use
  D=Get(nameIn);
  if (D!=NULL) {
    rprintf("  Error: DataSetList::Add: Data set %s already defined.\n",nameIn);
    return NULL;
  }
  switch (inType) {
    //case DOUBLE       : D = new dDataSet(); break;
    case DOUBLE       : D = new mapDataSet(); break;
    case STRING       : D = new stringDataSet(); break;
    case INT          : D = new intDataSet(); break;
    case UNKNOWN_DATA :
    default           :
      rprintf("  Error: DataSetList::Add: Unknown set type.\n");
      return NULL;
  }
  if (D==NULL) return NULL;
  // Set up dataset
  if ( D->Setup(nameIn,maxFrames) ) {
    rprintf("  Error setting up data set %s.\n",nameIn);
    delete D;
    return NULL;
  }

  DataList=(DataSet**) realloc(DataList,(Ndata+1) * sizeof(DataSet*));
  DataList[Ndata++]=D;
  //fprintf(stderr,"ADDED dataset %s\n",nameIn);
  return D;
}

/*
 * DataSetList::AddData()
 * Add data to a specific dataset in the list
 * Return 1 on error.
 */
int DataSetList::AddData(int frame, void *dataIn, int SetNumber) {
  if (SetNumber<0 || SetNumber>=Ndata) return 1;
  DataList[SetNumber]->Add(frame, dataIn);
  return 0;
} 

/*
 * DataSetList::Info()
 * Print information on all data sets in the list, as well as any datafiles
 * that will be written to.
 */
void DataSetList::Info() {
  int ds;

  mprintf("\nDATASETS:\n");
  if (Ndata==0)
    mprintf("  There are no data sets.");
  else if (Ndata==1)
    mprintf("  There is 1 data set: ");
  else
    mprintf("  There are %i data sets: ",Ndata);

  for (ds=0; ds<Ndata; ds++) {
    if (ds>0) mprintf(",");
    mprintf("%s",DataList[ds]->Name());
    //DataList[ds]->Info();
  }
  mprintf("\n");

  // DataFile Info
  //DFL.Info();
}

/*
 * DataSetList::Sync()
 * Call Sync for all datasets in the list
 */
void DataSetList::Sync() {
  int ds;

  // Sync datasets - does nothing if worldsize is 1
  for (ds=0; ds<Ndata; ds++) {
    if ( DataList[ds]->Sync() ) {
      rprintf( "Error syncing dataset %i\n",ds);
      //return;
    }
  }

}
