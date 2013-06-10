// DataFileList
#include "DataFileList.h"
#include "CpptrajStdio.h"
#ifdef MPI
# include "StringRoutines.h" // integerToString
#endif

// CONSTRUCTOR
DataFileList::DataFileList() : 
  debug_(0)
#ifdef MPI
  ,ensembleMode_(-1)
#endif
{}

// DESTRUCTOR
DataFileList::~DataFileList() {
  Clear();
}

// DataFileList::Clear()
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
  for (DFarray::iterator df = fileList_.begin(); df != fileList_.end(); ++df)
    (*df)->SetDebug( debug_ );
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
// TODO: Accept const ArgList so arguments are not reset?
DataFile* DataFileList::AddDataFile(std::string const& nameIn, ArgList& argIn) {
  // If no filename, no output desired
  if (nameIn.empty()) return 0;
  std::string name = nameIn;
# ifdef MPI
  if (ensembleMode_ != -1)
    // Ensemble mode, append rank to the output filename.
    name += ("." + integerToString(ensembleMode_));
# endif
  // Check if this filename already in use
  DataFile* Current = GetDataFile(name);
  // If no DataFile associated with name, create new datafile
  if (Current==0) {
    Current = new DataFile();
    if (Current->SetupDatafile(name, argIn, debug_)) {
      mprinterr("Error setting up DataFile %s\n",name.c_str());
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
/** Add given DataSet to the specified DataFile. If the DataFile does not
  * exist it will be created. Whenever a set is added to a data file
  * reset its writeFile status to true.
  */
DataFile* DataFileList::AddSetToFile(std::string const& nameIn, DataSet* dsetIn) {
  DataFile* DF = AddDataFile( nameIn );
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

// DataFileList::WriteAllDF()
/** Call write for all DataFiles in list for which writeFile is true. Once
  * a file has been written set writeFile to false; it can be reset to
  * true if new DataSets are added to it.
  */
void DataFileList::WriteAllDF() {
  for (DFarray::iterator df = fileList_.begin(); df != fileList_.end(); ++df) {
    mprintf("DBG: File %s writeFile=%i\n", (*df)->DataFilename().base(), (int)(*df)->DFLwrite());
    if ( (*df)->DFLwrite() ) {
      (*df)->Write();
      (*df)->SetDFLwrite( false );
    }
  }
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
