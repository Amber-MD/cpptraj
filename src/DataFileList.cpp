// DataFileList
#include <cstddef> // NULL
#include "DataFileList.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataFileList::DataFileList() :
  debug_(0)
{
//  fprintf(stderr,"DataFileList CONSTRUCTOR\n");
}

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
DataFile *DataFileList::GetDataFile(const char *nameIn) {
  if (nameIn==NULL) return NULL;
  int idx = FindName( nameIn );
  if (idx == -1) return NULL;
  return fileList_[idx];
}

// DataFileList::Add()
/** Add dataset to datafile in list with given file name. If the file does
  * not yet exist in the list create it. Return a pointer to the datafile
  * in the list.
  */
DataFile *DataFileList::Add(const char *nameIn, DataSet *dsetIn) {
  // If no filename, no output desired
  if (nameIn==NULL) return NULL;
  // If DataSet is NULL, dont add
  if (dsetIn==NULL) {
    mprintf("Error: Attempting to add non-existent dataset to file %s\n",nameIn);
    return NULL;
  }

  // Check if this filename already in use
  DataFile *Current = GetDataFile(nameIn);

  // If no DataFile associated with nameIn, create new datafile
  if (Current==NULL) {
    Current = new DataFile();
    if (Current->SetupDatafile(nameIn)) {
      mprinterr("Error setting up DataFile %s\n",nameIn);
      delete Current;
      return NULL;
    } 
    fileList_.push_back(Current);
    AddFilename( (char*)nameIn );
  }

  // Add the dataset to the current DataFile
  Current->AddSet(dsetIn);

  // Set debug level
  Current->SetDebug(debug_);

  // DEBUG
  //mprintf("** ADDED DATASET %s TO FILE %s\n",D->Name(),Current->filename);

  return Current;
}

// DataFileList::Info()
/** Print information on what datasets are going to what datafiles */
void DataFileList::Info() {
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
void DataFileList::ProcessDataFileArgs(DataSetList *masterDSL) {
  std::string df_cmd;
  char *name1 = NULL;
  char *name2 = NULL;
  DataFile *df;

  if (DF_Args_.empty()) return;
  mprintf("DATAFILE SETUP:\n");

  // Loop through all "datafile" arguments
  for (std::vector<ArgList>::iterator dataArg = DF_Args_.begin();
                                      dataArg != DF_Args_.end();
                                      dataArg++)
  {
    // Next string will be the argument passed to datafile
    df_cmd = (*dataArg).GetStringNext();
    if (df_cmd.empty()) {
      mprintf("Warning: datafile: No command given.\n");
      continue;
    }
    mprintf("  [%s]\n",(*dataArg).ArgLine());
    // Process command
    if ( df_cmd == "create" ) {
      // Usage: datafile create <filename> <dataset0> <dataset1> ...
      // Next string is datafile that command pertains to.
      name1 = (*dataArg).getNextString();
      if (name1==NULL) {
        mprinterr("Error: datafile create: No filename given.\n");
        mprinterr("Error: Usage: datafile create <filename> <set1> [<set2> ...]\n");
        continue;
      }
      df = GetDataFile(name1);
      if (df==NULL)
        mprintf("    Creating file %s\n",name1);
      else
        mprintf("    Adding sets to file %s\n",name1);
      while ( (name2=(*dataArg).getNextString())!=NULL ) {
        if ( Add(name1, masterDSL->Get(name2))==NULL ) {
          mprintf("Warning: Dataset %s does not exist in main dataset list, skipping.\n",name2);
        }
      }
    } else if ( df_cmd == "precision" ) {
      // Usage: datafile precision <filename> <dataset> [<width>] [<precision>]
      //        If width/precision not specified default to 12.4
      // Next string is datafile that command pertains to.
      name1 = (*dataArg).getNextString();
      if (name1==NULL) {
        mprinterr("Error: datafile precision: No filename given.\n");
        mprinterr("Error: Usage: datafile precision <filename> [<set>] <width> [<precision>]\n");
        continue;
      }
      df = GetDataFile(name1);
      if (df==NULL) {
        mprinterr("Error: datafile precision: DataFile %s does not exist.\n",name1);
        continue;
      }
      // This will break if dataset name starts with a digit...
      int width = (*dataArg).getNextInteger(12);
      int precision = (*dataArg).getNextInteger(4);
      name2 = (*dataArg).getNextString();
      df->SetPrecision(name2,width,precision);
    } else {
      // Argument is not a command but a datafile name that 
      // args will be passed to.
      df = GetDataFile( df_cmd.c_str() );
      if (df == NULL) {
        mprinterr("Error: datafile: File %s not found.\n", df_cmd.c_str());
        continue;
      }
      // Pass the rest of the arglist to the datafile.
      // FIXME: Exit on errors?
      df->ProcessArgs( *dataArg );
    }
  } // END loop over datafile args
}

// DataFileList::AddDatafileArg()
void DataFileList::AddDatafileArg(ArgList &argIn) {
  DF_Args_.push_back( argIn );
}

