// DataSetList
#include <cstdio> // sprintf
#include <cstdlib>
#include <cstring>
// This also includes basic DataSet class and dataType
#include "DataSetList.h"
#include "CpptrajStdio.h"
// Data types go here
#include "doubleDataSet.h"
#include "stringDataSet.h"
#include "intDataSet.h"
#include "mapDataSet.h"
#include "DataSet_float.h"

// CONSTRUCTOR
DataSetList::DataSetList() {
  //fprintf(stderr,"DSL Constructor\n");
  DataList=NULL;
  Ndata=0;
  maxFrames=0;
  debug=0;
  currentSet=0;
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

/* DataSetList::SetDebug()
 */
void DataSetList::SetDebug(int debugIn) {
  debug = debugIn;
  if (debug>0) mprintf("DataSetList Debug Level set to %i\n",debug);
}

/* DataSetList::SetMax()
 * Set the max number frames expected to be read in. Used to preallocate
 * data set sizes in the list.
 */
void DataSetList::SetMax(int expectedMax) {
  maxFrames = expectedMax;
  if (maxFrames<0) maxFrames=0;
}

/* DataSetList::Get()
 * Return pointer to DataSet with given name
 */
DataSet *DataSetList::Get(char *nameIn) {
  for (int i=0; i<Ndata; i++)
    if ( strcmp(DataList[i]->Name(),nameIn)==0 ) return DataList[i];
  return NULL;
}

/* DataSetList::Get()
 * Return pointer to DataSet with given idx
 */
DataSet *DataSetList::Get(int idxIn) {
  for (int i=0; i<Ndata; i++)
    if (DataList[i]->Idx() == idxIn) return DataList[i];
  return NULL;
}

/* DataSetList::AddMult()
 * Works like Add, except the dataset will be named 
 *   nameprefix_namesuffx
 * if nameprefix is not NULL, and
 *   namesuffix_XXXXX
 * otherwise.
 */
DataSet *DataSetList::AddMulti(dataType inType, char *prefix, const char *suffix) {
  char tempName[32];
  if (prefix==NULL)
    return this->Add(inType,NULL,suffix);
 // Sanity check - make sure name is not too big
 if (strlen(prefix)+strlen(suffix)+1 > 32) {
   mprinterr("Internal Error: DataSetList::AddMulti: size of %s+%s > 32\n",prefix,suffix);
   return NULL;
 } 
 sprintf(tempName,"%s_%s",prefix,suffix);
 return this->Add(inType,tempName,suffix);
}

/* DataSetList::Add()
 * Add a DataSet of specified type, set it up and return pointer to it. 
 * Name the dataset nameIn if specified, otherwise give it a default
 * name based on the given defaultName and dataset #. This routine
 * MUST be called with a default name.
 */ 
DataSet *DataSetList::Add(dataType inType, char *nameIn, const char *defaultName) {
  DataSet *D=NULL;
  char tempName[32];

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
    mprinterr("Error: DataSetList::Add: Data set %s already defined.\n",nameIn);
    return NULL;
  }
  switch (inType) {
    case DOUBLE       : D = new doubleDataSet(); break;
    case FLOAT        : D = new DataSet_float(); break;
    case STRING       : D = new stringDataSet(); break;
    case INT          : D = new intDataSet(); break;
    case MAP          : D = new mapDataSet(); break;
    case UNKNOWN_DATA :
    default           :
      mprinterr("Error: DataSetList::Add: Unknown set type.\n");
      return NULL;
  }
  if (D==NULL) return NULL;
  // Set up dataset
  if ( D->Setup(nameIn,maxFrames) ) {
    mprinterr("Error setting up data set %s.\n",nameIn);
    delete D;
    return NULL;
  }

  DataList=(DataSet**) realloc(DataList,(Ndata+1) * sizeof(DataSet*));
  DataList[Ndata++]=D;
  //fprintf(stderr,"ADDED dataset %s\n",nameIn);
  return D;
}

/* DataSetList::AddIdx()
 * Add a dataset to the list with given name and specific numeric index. The
 * intended use is for things like setting up data for a list of residues,
 * where the residues may not be sequential or start from 0. Since this 
 * routine is intended for use internally (i.e. the residue names and numbers
 * are generated inside other functions like DSSP and PerResRMSD) only
 * print warnings for higher debug levels. 
 */
DataSet *DataSetList::AddIdx(dataType inType, char *nameIn, int idxIn) {
  DataSet *D = NULL;

  // Check if dataset name is already in use
  D=Get(nameIn);
  if (D!=NULL) {
    if (debug>0) 
      mprintf("Warning: DataSetList::AddIdx: Data set %s already defined.\n",nameIn);
    return NULL;
  }
  // Check if dataset index already in use
  D=Get(idxIn);
  if (D!=NULL) {
    if (debug>0) 
      mprintf("Warning: DataSetList::AddIdx: Data set index %i already defined.\n",idxIn);
    return NULL;
  }

  // Allocate dataset type
  switch (inType) {
    case DOUBLE       : D = new doubleDataSet(); break;
    case FLOAT        : D = new DataSet_float(); break;
    case STRING       : D = new stringDataSet(); break;
    case INT          : D = new intDataSet(); break;
    case MAP          : D = new mapDataSet(); break;
    case UNKNOWN_DATA :
    default           :
      mprinterr("Error: DataSetList::AddIdx: Unknown set type.\n");
      return NULL;
  }
  if (D==NULL) return NULL;

  // Set up dataset
  if ( D->Setup(nameIn,maxFrames) ) {
    mprinterr("Error: DataSetList::AddIdx: setting up data set %s (%i).\n",nameIn,idxIn);
    delete D;
    return NULL;
  }
  D->SetIdx(idxIn);

  DataList=(DataSet**) realloc(DataList,(Ndata+1) * sizeof(DataSet*));
  DataList[Ndata++]=D;
  //fprintf(stderr,"ADDED dataset %s\n",nameIn);
  return D;
}

/* DataSetList::Begin()
 * Reset the set counter to 0.
 */
void DataSetList::Begin() {
  currentSet=0;
}

/* DataSetList::AddData()
 * Add data to the currentSet and increment the counter. Return 1 if at the 
 * last set. Should be used in conjunction with Begin. 
 */
int DataSetList::AddData(int frame, void *dataIn) {
  if (currentSet == Ndata) return 1;
  DataList[currentSet]->Add(frame, dataIn);
  currentSet++;
  return 0;
}

/* DataSetList::AddData()
 * Add data to a specific dataset in the list
 * Return 1 on error.
 */
int DataSetList::AddData(int frame, void *dataIn, int SetNumber) {
  if (SetNumber<0 || SetNumber>=Ndata) return 1;
  DataList[SetNumber]->Add(frame, dataIn);
  return 0;
}

/* DataSetList::AddDataToIdx
 */
/*
int DataSetList::AddDataToIdx(int frame, void *dataIn, int idxIn) {
  DataSet *D = this->Get(idxIn);
  if (D!=NULL) {
    D->Add(frame, dataIn);
    return 0;
  }
  return 1;
}*/

/* DataSetList::Info()
 * Print information on all data sets in the list, as well as any datafiles
 * that will be written to.
 */
void DataSetList::Info() {
  int ds;

  mprintf("\nDATASETS:\n");
  if (Ndata==0)
    mprintf("  There are no data sets set up for analysis.");
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

/* DataSetList::Sync()
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

