// DataSetList
#include <cstdio> // sprintf
#include <cstring>
//#include <stdexcept>
// This also includes basic DataSet class and dataType
#include "DataSetList.h"
#include "CpptrajStdio.h"
// Data types go here
#include "DataSet_double.h"
#include "DataSet_string.h"
#include "DataSet_integer.h"
#include "DataSet_XYZ.h"
#include "DataSet_float.h"
#include "DataSet_Matrix.h"

// CONSTRUCTOR
DataSetList::DataSetList() {
  //fprintf(stderr,"DSL Constructor\n");
  Ndata=0;
  maxFrames=0;
  debug=0;
  //currentSet=0;
}

// DESTRUCTOR
DataSetList::~DataSetList() {
  //fprintf(stderr,"DSL Destructor\n"); 
    for (int i=0; i<Ndata; i++)
      delete DataList[i];
}

// DataSetList::SetDebug()
void DataSetList::SetDebug(int debugIn) {
  debug = debugIn;
  if (debug>0) mprintf("DataSetList Debug Level set to %i\n",debug);
}

// DataSetList::SetMax()
/** Set the max number frames expected to be read in. Used to preallocate
  * data set sizes in the list.
  */
void DataSetList::SetMax(int expectedMax) {
  maxFrames = expectedMax;
  if (maxFrames<0) maxFrames=0;
}

/* DataSetList::SetPrecisionOfDatasets()
 * Set the width and precision for all datasets in the list.
 */
void DataSetList::SetPrecisionOfDatasets(int widthIn, int precisionIn) {
  for (int ds = 0; ds < Ndata; ds++) 
    DataList[ds]->SetPrecision(widthIn,precisionIn);
}

// DataSetList::operator[]()
/*DataSet &DataSetList::operator[](int ndataset) {
  if (ndataset < 0 || ndataset >= Ndata)
    throw std::out_of_range("DataSetList[]");
  return *(DataList[ndataset]);
}*/

// DataSetList::Get()
DataSet *DataSetList::Get(char *nameIn) {
  for (int i=0; i<Ndata; i++)
    if ( strcmp(DataList[i]->Name(),nameIn)==0 ) return DataList[i];
  return NULL;
}

// DataSetList::GetDataSetIdx()
DataSet *DataSetList::GetDataSetIdx(int idxIn) {
  for (int i=0; i<Ndata; i++)
    if (DataList[i]->Idx() == idxIn) return DataList[i];
  return NULL;
}

// DataSetList::GetDataSetN()
DataSet *DataSetList::GetDataSetN(int ndataset) {
  if (ndataset < 0 || ndataset >= Ndata) return NULL;
  return DataList[ndataset];
}

// DataSetList::AddMultiN()
/** Add a dataset to the list, basing the dataset name on prefix,
  * suffix, and a given number. If a dataset already exists in
  * the list NULL will be returned instead. The dataset will be named
  *   <prefix>_<suffix><Nin> 
  * if <prefix> is not NULL, 
  *   <suffix><Nin> 
  * if <prefix> is blank (i.e. ""), and 
  *   <suffix><Nin>_XXXXX
  * otherwise. The intended use is for things like setting up data for a 
  * list of residues, where the residues may not be sequential or start 
  * from 0. 
  * \param inType type of dataset to set up
  * \param prefix dataset name prefix
  * \param suffix dataset name suffix
  * \param Nin Number that can be used to uniquely identify the dataset.
  * \return pointer to successfully set up dataset.
  */
// NOTE: Instead of using Nin to identify, why not use dataset name?
DataSet *DataSetList::AddMultiN(dataType inType, const char *prefix, 
                                const char *suffix, int Nin) {
  char tempName[64];
  char tempSuffix[32];
  DataSet *tempDS = NULL;
  // Determine if dataset with idx Nin exists.
  tempDS = GetDataSetIdx( Nin );
  if (tempDS != NULL) return NULL; 
  sprintf(tempSuffix,"%s%i",suffix,Nin);
  if (prefix==NULL) {
    tempDS = this->Add(inType,NULL,tempSuffix);
  } else if (strcmp(prefix,"")==0) {
    tempDS = this->Add(inType,tempSuffix,tempSuffix);
  // Sanity check - make sure name is not too big
  } else {
    if (strlen(prefix)+strlen(tempSuffix)+1 > 64) {
      mprinterr("Internal Error: DataSetList::AddMultiN: size of %s+%s > 64\n",prefix,tempSuffix);
      return NULL;
    }
    sprintf(tempName,"%s_%s",prefix,tempSuffix);
    tempDS = this->Add(inType,tempName,tempSuffix);
  }
  if (tempDS!=NULL) tempDS->SetIdx( Nin );
  return tempDS;
}

// DataSetList::AddMulti()
/** Works like Add, except the dataset will be named 
  *   <prefix>_<suffix>
  * if <prefix> is not NULL, and
  *   <suffix>_XXXXX
  * otherwise.
  * \param inType type of dataset to set up
  * \param prefix dataset name prefix
  * \param suffix dataset name suffix
  * \return pointer to successfully set up dataset.
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

// DataSetList::checkName()
/** Given a potential dataset name, check if that name exists in the
  * DataSetList. If no name is given, generate a name based on
  * a default name and the dataset number.
  * \return An acceptable dataset name, or NULL if not acceptable.
  */
char *DataSetList::checkName(char *nameIn, const char *defaultName) {
  char *dsetName;
  size_t namesize;
  // Require all calls provide a default name
  if (defaultName==NULL) {
    mprinterr("Internal Error: DataSetList called without default name.\n");
    return NULL;
  }
  // If no name given make one up based on the given default 
  if (nameIn==NULL) {
    // Determine size of name + extension
    namesize = strlen( defaultName );
    size_t extsize = (size_t) DigitWidth( Ndata ); // # digits
    if (extsize < 5) extsize = 5;                  // Minimum size is 5 digits
    extsize += 2;                                  // + underscore + null
    namesize += extsize;
    dsetName = new char[ namesize ];
    sprintf(dsetName,"%s_%05i",defaultName,Ndata);
    //mprintf("NAME SIZE [%s] = %lu\n",dsetName,namesize);
  } else {
    namesize = strlen( nameIn ) + 1; // + null
    dsetName = new char[ namesize ];
    strcpy(dsetName, nameIn);
  }
  // Check if dataset name is already in use
  if (Get(dsetName)!=NULL) {
    mprinterr("Error: Data set %s already defined.\n",dsetName);
    delete[] dsetName;
    return NULL;
  }
  // Otherwise this dataset name is good to go
  return dsetName;
}

// DataSetList::AddMatrix()
/** Add a matrix dataset to the list. This requires a separate function
  * since the dimensions of the matrix must be set.
  */
DataSet *DataSetList::AddMatrix(char *nameIn, const char *defaultName,
                                int Nrows, int Mcols)
{
  DataSet *Dset=NULL;
  char *dsetName = checkName(nameIn, defaultName);
  if (dsetName==NULL) return NULL;
  Dset = new DataSet_Matrix(Nrows, Mcols);
  if (Dset==NULL) {
    delete[] dsetName;
    return NULL;
  }
  // Set up dataset
  if ( Dset->Setup(dsetName,maxFrames) ) {
    mprinterr("Error setting up matrix data set %s.\n",dsetName);
    delete Dset;
    delete[] dsetName;
    return NULL;
  }

  DataList.push_back(Dset);
  ++Ndata;
  //fprintf(stderr,"ADDED dataset %s\n",dsetName);
  delete[] dsetName;
  return Dset;
}

// DataSetList::Add()
/** Add a DataSet of specified type, set it up and return pointer to it. 
  * Name the dataset nameIn if specified, otherwise give it a default
  * name based on the given defaultName and dataset #. This routine
  * MUST be called with a default name.
  * \param inType type of dataset to set up
  * \param nameIn dataset name, can be NULL.
  * \param defaultName default name prefix for use if nameIn not specified.
  * \return pointer to successfully set-up dataset.
  */ 
DataSet *DataSetList::Add(dataType inType, char *nameIn, const char *defaultName) {
  DataSet *D=NULL;
  char *dsetName = checkName(nameIn, defaultName);
  if (dsetName==NULL) return NULL;

  switch (inType) {
    case DOUBLE       : D = new DataSet_double(); break;
    case FLOAT        : D = new DataSet_float(); break;
    case STRING       : D = new DataSet_string(); break;
    case INT          : D = new DataSet_integer(); break;
    case XYZ          : D = new DataSet_XYZ(); break;
    case UNKNOWN_DATA :
    default           :
      mprinterr("Error: DataSetList::Add: Unknown set type.\n");
      return NULL;
  }
  if (D==NULL) {
    delete[] dsetName;
    return NULL;
  }
  // Set up dataset
  if ( D->Setup(dsetName,maxFrames) ) {
    mprinterr("Error setting up data set %s.\n",dsetName);
    delete D;
    delete[] dsetName;
    return NULL;
  }

  DataList.push_back(D); 
  ++Ndata;
  //fprintf(stderr,"ADDED dataset %s\n",dsetName);
  delete[] dsetName;
  return D;
}

// DataSetList::AddData()
/** Add data to a specific dataset in the list
  * \param frame frame to add data to
  * \param dataIn pointer to data to add
  * \param SetNumber dataset to add data to
  * \return 0 on success, 1 on error.
  */
int DataSetList::AddData(int frame, void *dataIn, int SetNumber) {
  if (SetNumber<0 || SetNumber>=Ndata) return 1;
  DataList[SetNumber]->Add(frame, dataIn);
  return 0;
}

// DataSetList::Info()
/** Print information on all data sets in the list, as well as any datafiles
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

// DataSetList::Sync()
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

