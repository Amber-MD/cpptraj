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
    (*df)->RemoveDataSet( dsIn );
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
DataFile* DataFileList::GetDataFile(FileName const& nameIn) const {
  if (!nameIn.empty()) {
    for (DFarray::const_iterator df = fileList_.begin(); df != fileList_.end(); ++df)
      if (nameIn.Full() == (*df)->DataFilename().Full()) return *df;
  }
  return 0;
}

/** \return Index of CpptrajFile with specified file name if it exists in 
  *         the list, otherwise return -1. Must match full path.
  */
int DataFileList::GetCpptrajFileIdx(FileName const& nameIn) const {
  if (!nameIn.empty()) {
    for (int idx = 0; idx != (int)cfList_.size(); idx++)
      if (nameIn.Full() == cfList_[idx]->Filename().Full()) return idx;
  }
  return -1;
}

CpptrajFile* DataFileList::GetCpptrajFile(FileName const& nameIn) const {
  int idx = GetCpptrajFileIdx( nameIn );
  if (idx == -1) return 0;
  return cfList_[idx];
}

/** This version does not reset incoming arguments. */
DataFile* DataFileList::AddDataFile(FileName const& nameIn, DataFile::DataFormatType typeIn,
                                    ArgList const& argIn)
{
  ArgList args( argIn );
  return AddDataFile( nameIn, args, typeIn );
}

/** Create new DataFile, or return existing DataFile. */
DataFile* DataFileList::AddDataFile(FileName const& nameIn, ArgList& argIn,
                                    DataFile::DataFormatType typeIn)
{
  // If no filename, no output desired
  if (nameIn.empty()) return 0;
  // Check if filename in use by CpptrajFile.
  CpptrajFile* cf = GetCpptrajFile(nameIn);
  if (cf != 0) {
    mprinterr("Error: Data file name '%s' already in use by text output file '%s'.\n",
              nameIn.full(), cf->Filename().full());
    return 0;
  }
  // Check if this filename already in use
  DataFile* Current = GetDataFile(nameIn);
  // If no DataFile associated with name, create new DataFile
  if (Current==0) {
    Current = new DataFile();
    if (Current->SetupDatafile(nameIn, argIn, typeIn, debug_)) {
      mprinterr("Error: Setting up data file %s\n", nameIn.full());
      delete Current;
      return 0;
    }
    fileList_.push_back(Current);
  } else {
    // Set debug level
    Current->SetDebug(debug_);
    // If a type was specified, make sure it matches.
    if (typeIn != DataFile::UNKNOWN_DATA && typeIn != Current->Type()) {
      mprinterr("Error: '%s' is type %s but has been requested as type %s.\n",
                Current->DataFilename().full(), Current->FormatString(),
                DataFile::FormatString( typeIn ));
      return 0;
    }
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

DataFile* DataFileList::AddDataFile(FileName const& nameIn, ArgList& argIn) {
  return AddDataFile( nameIn, argIn, DataFile::UNKNOWN_DATA );
}

// DataFileList::AddDataFile()
DataFile* DataFileList::AddDataFile(FileName const& nameIn) {
  ArgList empty;
  return AddDataFile( nameIn, empty, DataFile::UNKNOWN_DATA );
}

// DataFileList::AddCpptrajFile()
/** File type is text, stdout not allowed. */
CpptrajFile* DataFileList::AddCpptrajFile(FileName const& nameIn, std::string const& descrip)
{ return AddCpptrajFile(nameIn, descrip, TEXT, false); }

// DataFileList::AddCpptrajFile()
/** Stdout is not allowed. */
CpptrajFile* DataFileList::AddCpptrajFile(FileName const& nameIn, 
                                          std::string const& descrip, CFtype typeIn)
{ return AddCpptrajFile(nameIn, descrip, typeIn, false); }

/** Create new CpptrajFile of the given type, or return existing CpptrajFile.
  * STDOUT will be used if name is empty and STDOUT is allowed.
  */
// TODO: Accept const ArgList so arguments are not reset?
CpptrajFile* DataFileList::AddCpptrajFile(FileName const& nameIn, 
                                          std::string const& descrip,
                                          CFtype typeIn, bool allowStdout)
{
  // If no filename and stdout not allowed, no output desired.
  if (nameIn.empty() && !allowStdout) return 0;
  FileName name;
  CpptrajFile* Current = 0;
  int currentIdx = -1;
  if (!nameIn.empty()) {
#   ifdef MPI
    // FIXME: Unlike DataFiles, CpptrajFiles are opened immediately so append
    //        worldrank to filename now. This will have to change if MPI does
    //        not necessarily mean ensemble mode in the future.
    name.SetFileName( AppendNumber(nameIn.Full(), worldrank) );
#   else
    name = nameIn;
#   endif
    // Check if filename in use by DataFile.
    DataFile* df = GetDataFile(name);
    if (df != 0) {
      mprinterr("Error: Text output file name '%s' already in use by data file '%s'.\n",
                nameIn.full(), df->DataFilename().full());
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
      mprinterr("Error: Setting up text output file %s\n", name.full());
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
    mprintf("\nDATAFILES (%zu total):\n", fileList_.size() + cfList_.size());
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
      (*df)->WriteDataOut();
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
