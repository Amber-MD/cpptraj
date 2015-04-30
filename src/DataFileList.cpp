// DataFileList
#include "DataFileList.h"
#include "CpptrajStdio.h"
#include "MpiRoutines.h" // parallel_barrier
#include "PDBfile.h"
#ifdef MPI
# include "StringRoutines.h" // integerToString
#endif
#ifdef TIMER
# include "Timer.h"
#endif

// CONSTRUCTOR
DataFileList::DataFileList() : debug_(0) {}

// DESTRUCTOR
DataFileList::~DataFileList() { Clear(); }

// DataFileList::Clear()
void DataFileList::Clear() {
  for (DFarray::iterator it = fileList_.begin(); it != fileList_.end(); ++it)
    delete *it;
  fileList_.clear();
  for (CFarray::iterator it = cfList_.begin(); it != cfList_.end(); ++it) {
    (*it)->CloseFile();
    delete *it;
  }
  cfList_.clear();
  cfData_.clear();
}

// DataFileList::RemoveDataFile()
DataFile* DataFileList::RemoveDataFile( DataFile* dfIn ) {
  for (DFarray::iterator it = fileList_.begin(); it != fileList_.end(); ++it) {
    if ( dfIn == *it ) {
      delete *it;
      return (DataFile*)0;
    }
  }
  return dfIn;
}

// DataFileList::RemoveDataSet()
/** Remove given DataSet from any DataFiles in list. */
void DataFileList::RemoveDataSet( DataSet* dsIn ) {
  for (DFarray::iterator df = fileList_.begin(); df != fileList_.end(); ++df)
    (*df)->RemoveSet( dsIn );
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

#ifdef MPI
void DataFileList::MakeDataFilesEnsemble(int memberIn) {
  for (DFarray::const_iterator df = fileList_.begin(); df != fileList_.end(); ++df)
    (*df)->SetMember( memberIn );
}
#endif

// DataFileList::GetDataFile()
/** \return DataFile specified by given file name if it exists in the list,
  *         otherwise return 0. Must match full path.
  */
DataFile* DataFileList::GetDataFile(std::string const& nameIn) const {
  if (!nameIn.empty()) {
    for (DFarray::const_iterator df = fileList_.begin(); df != fileList_.end(); ++df)
      if (nameIn == (*df)->DataFilename().Full()) return *df;
  }
  return 0;
}

/** \return Index of CpptrajFile with specified file name if it exists in 
  *         the list, otherwise return -1. Must match full path.
  */
int DataFileList::GetCpptrajFileIdx(std::string const& nameIn) const {
  if (!nameIn.empty()) {
    for (int idx = 0; idx != (int)cfList_.size(); idx++)
      if (nameIn == cfList_[idx]->Filename().Full()) return idx;
  }
  return -1;
}

CpptrajFile* DataFileList::GetCpptrajFile(std::string const& nameIn) const {
  int idx = GetCpptrajFileIdx( nameIn );
  if (idx == -1) return 0;
  return cfList_[idx];
}

/** Create new DataFile, or return existing DataFile. */
// TODO: Accept const ArgList so arguments are not reset?
DataFile* DataFileList::AddDataFile(std::string const& nameIn, ArgList& argIn) {
  // If no filename, no output desired
  if (nameIn.empty()) return 0;
  std::string name = nameIn;
  // Check if filename in use by CpptrajFile.
  CpptrajFile* cf = GetCpptrajFile(name);
  if (cf != 0) {
    mprinterr("Error: Data file name '%s' already in use by text output file '%s'.\n",
              nameIn.c_str(), cf->Filename().full());
    return 0;
  }
  // Check if this filename already in use
  DataFile* Current = GetDataFile(name);
  // If no DataFile associated with name, create new DataFile
  if (Current==0) {
    Current = new DataFile();
    if (Current->SetupDatafile(name, argIn, debug_)) {
      mprinterr("Error: Setting up data file %s\n",name.c_str());
      delete Current;
      return 0;
    }
    fileList_.push_back(Current);
  } else {
    // Set debug level
    Current->SetDebug(debug_);
    // Check for keywords that do not match file type
    DataFile::DataFormatType kType = DataFile::GetFormatFromArg( argIn );
    if (kType != DataFile::UNKNOWN_DATA && kType != Current->Type())
      mprintf("Warning: %s is type %s but type %s keyword specified; ignoring keyword.\n",
              Current->DataFilename().full(), Current->FormatString(),
              DataFile::FormatString( kType ));
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

// DataFileList::AddCpptrajFile()
/** File type is text, stdout not allowed. */
CpptrajFile* DataFileList::AddCpptrajFile(std::string const& nameIn, std::string const& descrip)
{ return AddCpptrajFile(nameIn, descrip, TEXT, false); }

// DataFileList::AddCpptrajFile()
/** Stdout is not allowed. */
CpptrajFile* DataFileList::AddCpptrajFile(std::string const& nameIn, 
                                          std::string const& descrip, CFtype typeIn)
{ return AddCpptrajFile(nameIn, descrip, typeIn, false); }

/** Create new CpptrajFile of the given type, or return existing CpptrajFile.
  * STDOUT will be used if name is empty and STDOUT is allowed.
  */
// TODO: Accept const ArgList so arguments are not reset?
CpptrajFile* DataFileList::AddCpptrajFile(std::string const& nameIn, 
                                          std::string const& descrip,
                                          CFtype typeIn, bool allowStdout)
{
  // If no filename and stdout not allowed, no output desired.
  if (nameIn.empty() && !allowStdout) return 0;
  std::string name("");
  CpptrajFile* Current = 0;
  int currentIdx = -1;
  if (!nameIn.empty()) {
    name = nameIn;
#   ifdef MPI
    // FIXME: Unlike DataFiles, CpptrajFiles are opened immediately so append
    //        worldrank to filename now. This will have to change if MPI does
    //        not necessarily mean ensemble mode in the future.
    name.append("." + integerToString(worldrank));
#   endif
    // Check if filename in use by DataFile.
    DataFile* df = GetDataFile(name);
    if (df != 0) {
      mprinterr("Error: Text output file name '%s' already in use by data file '%s'.\n",
                nameIn.c_str(), df->DataFilename().full());
      return 0;
    }
    // Check if this filename already in use
    currentIdx = GetCpptrajFileIdx( name );
    if (currentIdx != -1) Current = cfList_[currentIdx];
  }
  // If no CpptrajFile associated with name, create new CpptrajFile
  if (Current==0) {
    switch (typeIn) {
      case TEXT: Current = new CpptrajFile(); break;
      case PDB:  Current = (CpptrajFile*)(new PDBfile()); break;
    }
    Current->SetDebug(debug_);
    // Set up file for writing. 
    //if (Current->SetupWrite( name, debug_ ))
    if (Current->OpenWrite( name ))
    {
      mprinterr("Error: Setting up text output file %s\n",name.c_str());
      delete Current;
      return 0;
    }
    cfList_.push_back( Current );
    cfData_.push_back( CFstruct(descrip, typeIn) );
  } else {
    // If Current type does not match typeIn do not allow.
    if (typeIn != cfData_[currentIdx].Type()) {
      mprinterr("Error: Cannot change type of text output for '%s'.\n", Current->Filename().full());
      return 0;
    }
    Current->SetDebug(debug_);
    // Update description
    if (!descrip.empty())
      cfData_[currentIdx].UpdateDescrip( descrip );
  }
  return Current;
}

// DataFileList::List()
/** Print information on what datasets are going to what datafiles */
void DataFileList::List() const {
  parallel_barrier();
  if (!fileList_.empty() || !cfList_.empty()) {
    mprintf("\nDATAFILES:\n");
    if (!fileList_.empty()) {
      for (DFarray::const_iterator it = fileList_.begin(); it != fileList_.end(); ++it)
        rprintf("  %s (%s): %s\n",(*it)->DataFilename().base(), (*it)->FormatString(),
                (*it)->DataSetNames().c_str());
    }
    if (!cfList_.empty()) {
      for (unsigned int idx = 0; idx != cfList_.size(); idx++)
        rprintf("  %s (%s)\n", cfList_[idx]->Filename().base(), cfData_[idx].descrip());
    }
  }
}

// DataFileList::WriteAllDF()
/** Call write for all DataFiles in list for which writeFile is true. Once
  * a file has been written set writeFile to false; it can be reset to
  * true if new DataSets are added to it.
  */
void DataFileList::WriteAllDF() {
  if (fileList_.empty()) return;
# ifdef TIMER
  Timer datafile_time;
  datafile_time.Start();
# endif
  for (DFarray::iterator df = fileList_.begin(); df != fileList_.end(); ++df) {
    if ( (*df)->DFLwrite() ) {
      (*df)->WriteData();
      (*df)->SetDFLwrite( false );
    }
  }
# ifdef TIMER
  datafile_time.Stop();
  mprintf("TIME: Write of all data files took %.4f seconds.\n", datafile_time.Total());
# endif
}

/** Reset writeFile status for all files in list to true. */
void DataFileList::ResetWriteStatus() {
  for (DFarray::iterator df = fileList_.begin(); df != fileList_.end(); ++df)
    (*df)->SetDFLwrite( true );
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
  int err = df->ProcessArgs( dataArg );
  if (err != 0 || dataArg.CheckForMoreArgs()) return 1;
  return 0;
}
