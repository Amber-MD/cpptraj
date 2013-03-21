// DataFileList
#include "DataFileList.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataFileList::DataFileList() : debug_(0) {}

// DESTRUCTOR
DataFileList::~DataFileList() {
  Clear();
}

void DataFileList::Clear() {
  for (DFarray::iterator it = fileList_.begin(); it != fileList_.end(); it++)
    delete *it;
  fileList_.clear();
  FileList::Clear();
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
DataFile* DataFileList::GetDataFile(std::string const& nameIn) const {
  if (nameIn.empty()) return 0;
  int idx = FindName( nameIn );
  if (idx == -1) return 0;
  return fileList_[idx];
}

/** Create new DataFile, or return existing DataFile. */
DataFile* DataFileList::AddDataFile(std::string const& nameIn, ArgList& argIn) {
  // If no filename, no output desired
  if (nameIn.empty()) return 0;
  // Check if this filename already in use
  DataFile* Current = GetDataFile(nameIn);
  // If no DataFile associated with nameIn, create new datafile
  if (Current==0) {
    Current = new DataFile();
    if (Current->SetupDatafile(nameIn, argIn, debug_)) {
      mprinterr("Error setting up DataFile %s\n",nameIn.c_str());
      delete Current;
      return 0;
    }
    fileList_.push_back(Current);
    AddFilename( Current->DataFilename() );
  } else {
    // Set debug level
    Current->SetDebug(debug_);
    // Process Arguments
    if (!argIn.empty())
      Current->ProcessArgs( argIn );
  }
  return Current;
}

// DataFileList::AddDataFile()
DataFile* DataFileList::AddDataFile(std::string const& nameIn) {
  ArgList empty;
  return AddDataFile( nameIn, empty );
}

// DataFileList::AddSetToFile()
DataFile* DataFileList::AddSetToFile(std::string const& nameIn, DataSet* dsetIn) {
  ArgList empty;
  DataFile* DF = AddDataFile( nameIn, empty );
  if (DF == 0) return 0;
  DF->AddSet( dsetIn );
  return DF;
}

// DataFileList::List()
/** Print information on what datasets are going to what datafiles */
void DataFileList::List() const {
  if (fileList_.empty()) {
    //mprintf("NO DATASETS WILL BE OUTPUT\n");
    return;
  }

  mprintf("DATAFILE OUTPUT:\n");
  for (DFarray::const_iterator it = fileList_.begin(); it != fileList_.end(); it++) {
    mprintf("  %s: ",(*it)->DataFilename().base());
    (*it)->DataSetNames();
    mprintf("\n");
  }
}

// DataFileList::Write()
/** Call write for all datafiles in list. Only master should call this.
  */
void DataFileList::Write() {
  //if (worldrank!=0) return; 
  for (DFarray::iterator it = fileList_.begin(); it != fileList_.end(); it++)
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
  if (df == 0) {
    mprinterr("Error: datafile: File %s not found.\n", df_cmd.c_str());
    return 1;
  }
  // Process command
  return df->ProcessArgs( dataArg );
}
