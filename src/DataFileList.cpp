// DataFileList
#include <cstddef> // NULL
#include "DataFileList.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataFileList::DataFileList() :
  debug_(0)
{ }

// DESTRUCTOR
DataFileList::~DataFileList() {
//  fprintf(stderr,"DataFileList DESTRUCTOR\n");
  for (df_iterator it = fileList_.begin(); it != fileList_.end(); it++)
    delete *it;
}

// DataFileList::SetDebug()
/** Set debug level for DataFileList and all datafiles in it. */
void DataFileList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_>0)
    mprintf("DataFileList DEBUG LEVEL SET TO %i\n",debug_);
}

// DataFileList::GetDataFile()
/** Return DataFile specified by given file name if it exists in the list,
  * otherwise return null.
  */
DataFile* DataFileList::GetDataFile(std::string const& nameIn) {
  if (nameIn.empty()) return NULL;
  int idx = FindName( nameIn );
  if (idx == -1) return NULL;
  return fileList_[idx];
}

// DataFileList::Add()
DataFile* DataFileList::Add(const char* nameIn, DataSet* dsetIn) {
  if (nameIn == NULL) return NULL;
  return AddSetToFile( std::string(nameIn), dsetIn );
}

// DataFileList::AddSetToFile()
/** Add dataset to datafile in list with given file name. If the file does
  * not yet exist in the list create it. Return a pointer to the datafile
  * in the list.
  */
DataFile* DataFileList::AddSetToFile(std::string const& nameIn, DataSet* dsetIn) {
  // If no filename, no output desired
  if (nameIn.empty()) return NULL;
  // If DataSet is NULL, dont add
  if (dsetIn==NULL) {
    mprintf("Error: Attempting to add non-existent dataset to file %s\n",nameIn.c_str());
    return NULL;
  }

  // Check if this filename already in use
  DataFile* Current = GetDataFile(nameIn);

  // If no DataFile associated with nameIn, create new datafile
  if (Current==NULL) {
    Current = new DataFile();
    if (Current->SetupDatafile(nameIn)) {
      mprinterr("Error setting up DataFile %s\n",nameIn.c_str());
      delete Current;
      return NULL;
    } 
    fileList_.push_back(Current);
    AddFilename( nameIn );
  }

  // Add the dataset to the current DataFile
  Current->AddSet(dsetIn);

  // Set debug level
  Current->SetDebug(debug_);

  // DEBUG
  //mprintf("** ADDED DATASET %s TO FILE %s\n",D->Name(),Current->filename);

  return Current;
}

// DataFileList::List()
/** Print information on what datasets are going to what datafiles */
void DataFileList::List() {
  if (fileList_.empty()) {
    //mprintf("NO DATASETS WILL BE OUTPUT\n");
    return;
  }

  mprintf("DATAFILE OUTPUT:\n");
  for (df_iterator it = fileList_.begin(); it != fileList_.end(); it++) {
    mprintf("  %s: ",(*it)->Filename());
    (*it)->DataSetNames();
    mprintf("\n");
  }
}

// DataFileList::Write()
/** Call write for all datafiles in list. Only master should call this.
  */
void DataFileList::Write() {
  //if (worldrank!=0) return; 
  for (df_iterator it = fileList_.begin(); it != fileList_.end(); it++)
    (*it)->Write();
}

// DataFileList::ProcessDataFileArgs()
/** Process command relating to data files. */
int DataFileList::ProcessDataFileArgs(ArgList& dataArg) {
  // Next string is DataFile name that command will be passed to.
  std::string df_cmd = dataArg.GetStringNext();
  if (df_cmd.empty()) {
    mprintf("Warning: datafile: No filename given.\n");
    return 0;
  }
  // Check for deprecated commands
  if (df_cmd == "create" || df_cmd == "precision") 
    mprintf("Warning: 'datafile %s' is deprecated; use %s instead.\n", 
            df_cmd.c_str(), df_cmd.c_str());
  //mprintf("  [%s]\n",(*dataArg).ArgLine());
  DataFile* df = GetDataFile( df_cmd.c_str() );
  if (df == NULL) {
    mprinterr("Error: datafile: File %s not found.\n", df_cmd.c_str());
    return 1;
  }
  // Process command
  return df->ProcessArgs( dataArg );
}

