#include <cstdio> // sscanf
#include <algorithm> // sort
#include "DataIO_RemLog.h"
#include "CpptrajStdio.h" 
#include "ProgressBar.h"
#include "StringRoutines.h" // fileExists

// CONSTRUCTOR
DataIO_RemLog::DataIO_RemLog() : searchForLogs_(true)
{
  SetValid( DataSet::REMLOG );
}

// NOTE: Must match LogType
const char* DataIO_RemLog::LogDescription[] = {
  "Unknown", "Temperature", "Hamiltonian", "MultipleDim", "RXSGLD", "pH"
};

/// \return true if char pointer is null.
static inline bool IsNullPtr( const char* ptr ) { 
  if (ptr == 0) {
    mprinterr("Error: Could not read file.\n");
    return true;
  }
  return false;
}

// DataIO_RemLog::ReadRemlogHeader()
/** Read information from given replica log file header.
  * \return Number of exchanges expected in log file.
  * \param buffer Log file, should be opened.
  * \param type Will be set to type of log file.
  * \param Ndims For MREMD, expected # of dims from remd.dim file.
  */
int DataIO_RemLog::ReadRemlogHeader(BufferedLine& buffer, LogType& type, unsigned int Ndims)
const {
  int numexchg = -1;
  type = UNKNOWN;
  // Read the first line. Should be '# Replica Exchange log file'
  std::string line = buffer.GetLine();
  if (line.compare(0, 27, "# Replica Exchange log file") != 0) {
    mprinterr("Error: Expected '# Replica Exchange log file', got:\n%s\n", line.c_str());
    return -1;
  }
  // Read past metadata. Save expected number of exchanges.
  while (line[0] == '#') {
    line = buffer.GetLine();
    if (line.empty()) {
      mprinterr("Error: No exchanges in rem log.\n");
      return -1;
    }
    ArgList columns( line );
    // Each line should have at least 2 arguments
    if (columns.Nargs() > 1) {
      if (columns[1] == "exchange")
        break;
      else
      {
        //if (debug_ > 0) mprintf("\t%s", line.c_str());
        if (columns[1] == "numexchg")
          numexchg = columns.getNextInteger(-1);
        else if (columns[1] == "Dimension") {
          type = MREMD;
          int ndim = columns.getKeyInt("of", 0);
          if (ndim != (int)Ndims) {
            mprinterr("Error: # of dimensions in rem log %i != dimensions in remd.dim file (%u).\n",
                      ndim, Ndims);
            return -1;
          }
        } else if (columns[1] == "RXSGLD")
          type = RXSGLD;
      }
    }
    if (type == UNKNOWN && columns.hasKey("Rep#,")) {
      if (columns[2] == "Neibr#,") type = HREMD;
      else if (columns[2] == "Velocity") type = TREMD;
      else if (columns[2] == "N_prot," ) type = PHREMD;
    }
  }
  if (numexchg < 1) {
    mprinterr("Error: Invalid number of exchanges (%i) in rem log.\n");
    return -1;
  }
  if (type == UNKNOWN) {
    mprinterr("Error: Unknown replica log type.\n");
    return -1;
  }
  return numexchg;
}

// -----------------------------------------------------------------------------
// DataIO_RemLog::ReadRemdDimFile()
// TODO: Handle cases where groups are not in order.
/** Read dimension info from replica dimension file. Set up Groups and dims.
  * \return 0 if setup OK, 1 if error.
  * \param rd_name replica dimension file name.
  * \param GroupDims Output array that will contain group layout.
  */
int DataIO_RemLog::ReadRemdDimFile(FileName const& rd_name,
                                   DataSet_RemLog::GdimArray& GroupDims,
                                   ReplicaDimArray& DimTypes)
{
  int n_mremd_replicas = 0;
  typedef std::map<int,DataSet_RemLog::GroupArray> GroupMapType;
  typedef std::pair<GroupMapType::iterator,bool> GroupMapRet;
  typedef std::pair<int,DataSet_RemLog::GroupArray> GroupMapElt;
  BufferedLine rd_file;
  if (rd_file.OpenFileRead( rd_name )) {
    mprinterr("Error: Could not read remd dim file '%s'\n", rd_name.full());
    return 0;
  }
  const char* separators = " =,()";
  // Read dimension title
  const char* ptr = rd_file.Line();
  if (IsNullPtr( ptr )) return 0;
  // ptr Should end with a newline
  mprintf("\tReplica dimension file '%s' title: %s\n", rd_name.base(), ptr);
  // Read each &multirem section
  GroupDims.clear();
  DimTypes.clear();
  ArgList rd_arg;
  while (ptr != 0) {
    rd_arg.SetList( std::string(ptr), separators );
    if ( rd_arg[0] == "&multirem" ) {
      GroupMapType GroupMap;
      std::string desc;
      int n_replicas = 0;
      ReplicaDimArray::RemDimType exch_type = ReplicaDimArray::UNKNOWN;
      while (ptr != 0) {
        rd_arg.SetList( std::string(ptr), separators );
        if (rd_arg.CommandIs("&end") || rd_arg.CommandIs("/")) break;
        rd_arg.MarkArg(0);
        if ( rd_arg.CommandIs("exch_type") ) {
          if ( rd_arg.hasKey("TEMP") || rd_arg.hasKey("TEMPERATURE") )
            exch_type = ReplicaDimArray::TEMPERATURE;
          else if ( rd_arg.hasKey("HAMILTONIAN") || rd_arg.hasKey("HREMD") )
            exch_type = ReplicaDimArray::HAMILTONIAN;
          else {
            mprinterr("Error: Unrecognized exch_type: %s\n", rd_arg.ArgLine());
            return 0;
          }
        } else if ( rd_arg.CommandIs("group") ) {
          int group_num = rd_arg.getNextInteger(-1);
          if (group_num < 1) {
            mprinterr("Error: Invalid group number: %i\n", group_num);
            return 0;
          }
          //mprintf("\t\tGroup %i\n", group_num);
          std::vector<int> indices;
          int group_index = rd_arg.getNextInteger(-1);
          while (group_index != -1) {
            indices.push_back( group_index );
            n_replicas++;
            group_index = rd_arg.getNextInteger(-1);
          }
          // Set up partner array for this group
          DataSet_RemLog::GroupArray group;
          for (int i = 0; i < (int)indices.size(); i++) {
            int l_idx = i - 1;
            if (l_idx < 0) l_idx = (int)indices.size() - 1;
            int r_idx = i + 1;
            if (r_idx == (int)indices.size()) r_idx = 0;
            group.push_back( DataSet_RemLog::
                             GroupReplica(indices[l_idx], indices[i], indices[r_idx]) );
            //mprintf("\t\t\t%i - %i - %i\n", group.back().L_partner(),
            //        group.back().Me(), group.back().R_partner());
          }
          GroupMapRet ret = GroupMap.insert( GroupMapElt(group_num, group) );
          if (ret.second == false) {
            mprinterr("Error: Duplicate group # detected (%i)\n", group_num);
            return 0;
          }
        } else if ( rd_arg.CommandIs("desc") ) {
          desc = rd_arg.GetStringNext(); 
        }
        ptr = rd_file.Line();
      }
      // Place sorted groups for dim into GroupDimType
      DataSet_RemLog::GroupDimType Groups;
      for (GroupMapType::const_iterator it = GroupMap.begin(); it != GroupMap.end(); ++it)
        Groups.push_back( it->second );
      mprintf("\tDimension %zu: type '%s', description '%s', groups=%zu, replicas=%i\n", 
              GroupDims.size() + 1, ReplicaDimArray::dimType(exch_type), desc.c_str(), 
              Groups.size(), n_replicas);
      if (n_mremd_replicas == 0)
        n_mremd_replicas = n_replicas;
      else if (n_replicas != n_mremd_replicas) {
        mprinterr("Error: Number of MREMD replicas in dimension (%i) != number of\n"
                  "Error: MREMD replicas in first dimension (%i)\n", n_replicas, n_mremd_replicas);
        return 0;
      }
      GroupDims.push_back( Groups );
      DimTypes.AddRemdDimension( exch_type );
    }
    ptr = rd_file.Line();
  }
  if (GroupDims.empty()) {
    mprinterr("Error: No replica dimensions found.\n");
    return 0;
  }
  return n_mremd_replicas;
}

// -----------------------------------------------------------------------------
// DataIO_RemLog::SetupTemperatureMap()
/** buffer should be positioned at the first exchange. */
DataIO_RemLog::TmapType 
  DataIO_RemLog::SetupTemperatureMap(BufferedLine& buffer,
                                     std::vector<int>& CrdIdxs) const
{
  TmapType TemperatureMap;
  std::vector<TlogType> tList; // Hold temps and associated coord idxs
  TlogType tlog;
  CrdIdxs.clear();
  const char* ptr = buffer.Line();
  while (ptr != 0 && ptr[0] != '#') {
    // For temperature remlog create temperature map. 
    //mprintf("DEBUG: Temp0= %s", ptr+32);
    if ( sscanf(ptr, "%2i%*10f%*10f%*10f%10lf", &tlog.crdidx, &tlog.t0) != 2 ) {
      mprinterr("Error: could not read temperature from T-REMD log.\n"
                "Error: Line: %s", ptr);
      return TemperatureMap;
    }
    tList.push_back( tlog );
    ptr = buffer.Line();
  }
  // Sort temperatures and associated coord indices
  std::sort( tList.begin(), tList.end(), TlogType_cmp() );
  // Place sorted temperatures into map starting from replica index 1. Check
  // for duplicate temperatures. Also store the sorted coordinate indices.
  int repidx = 1;
  for (std::vector<TlogType>::const_iterator it = tList.begin();
                                             it != tList.end(); ++it)
  {
    mprintf("\t\tReplica %i => %f (crdidx= %i)\n", repidx, it->t0, it->crdidx); 
    if (it != tList.begin()) {
      if ( it->t0 == (it-1)->t0 ) {
        mprinterr("Error: duplicate temperature %.2f detected in T-REMD remlog\n", it->t0);
        TemperatureMap.clear();
        return TemperatureMap;
      }
    }
    TemperatureMap.insert(std::pair<double,int>(it->t0, repidx++));
    CrdIdxs.push_back( it->crdidx );
  }

  return TemperatureMap;
}

// DataIO_RemLog::Setup_pH_Map()
/** buffer should be positioned at the first exchange. */
DataIO_RemLog::TmapType 
  DataIO_RemLog::Setup_pH_Map(BufferedLine& buffer, std::vector<int>& CrdIdxs) const
{
  TmapType pH_Map;
  std::vector<TlogType> pList; // Hold temps and associated coord idxs
  TlogType plog;
  CrdIdxs.clear();
  const char* ptr = buffer.Line();
  while (ptr != 0 && ptr[0] != '#') {
    // For pH remlog create pH map. 
    // '(i6,x,i7,x,2f7.3,x,f8.4)'
    if ( sscanf(ptr, "%6i%*8i%*c%7lf", &plog.crdidx, &plog.t0) != 2 ) {
      mprinterr("Error: could not read pH from pH-REMD log.\n"
                "Error: Line: %s", ptr);
      return pH_Map;
    }
    //mprintf("DEBUG: pH= %g\n", plog.t0);
    pList.push_back( plog );
    ptr = buffer.Line();
  }
  // Sort pHs and associated coord indices
  std::sort( pList.begin(), pList.end(), TlogType_cmp() );
  // Place sorted pHs into map starting from replica index 1. Check
  // for duplicate pHs. Also store the sorted coordinate indices.
  int repidx = 1;
  for (std::vector<TlogType>::const_iterator it = pList.begin();
                                             it != pList.end(); ++it)
  {
    mprintf("\t\tReplica %i => %f (crdidx= %i)\n", repidx, it->t0, it->crdidx); 
    if (it != pList.begin()) {
      if ( it->t0 == (it-1)->t0 ) {
        mprinterr("Error: duplicate pH %.2f detected in pH-REMD remlog\n", it->t0);
        pH_Map.clear();
        return pH_Map;
      }
    }
    pH_Map.insert(std::pair<double,int>(it->t0, repidx++));
    CrdIdxs.push_back( it->crdidx );
  }

  return pH_Map;
}

// DataIO_RemLog::CountHamiltonianReps()
/** buffer should be positioned at the first exchange. */
int DataIO_RemLog::CountHamiltonianReps(BufferedLine& buffer) const {
  int n_replicas = 0;
  const char* ptr = buffer.Line();
  while (ptr != 0 && ptr[0] != '#') {
    ptr = buffer.Line();
    ++n_replicas;
  }
  return n_replicas;
}

// -----------------------------------------------------------------------------
/** Open replica logs for all dimensions. For MREMD (# dims > 1), check that
  * the number of exchanges reported in each file is the same and that each
  * is actually a MREMD log. Will leave all remlogs positioned at the first
  * exchange.
  * \return Expected number of exchanges.
  * \param buffer Array of replica log files (unopened), 1 per dimension.
  * \param dimLogs Array of replica log file names, 1 per dimension.
  * \param expectedType expected type of log(s).
  */
int DataIO_RemLog::OpenMremdDims(std::vector<BufferedLine>& buffer, Sarray const& dimLogs,
                                 LogType expectedType)
{
  // Sanity check
  if (buffer.size() != dimLogs.size()) {
    mprinterr("Internal Error: File buffer array size %zu != # MREMD logs %zu.\n",
              buffer.size(), dimLogs.size());
    return 1;
  }
  int total_exchanges = -1;
  // Open remlogs for each dimension as buffered file.
  for (unsigned int dim = 0; dim != buffer.size(); dim++) {
    LogType log_type;
    buffer[dim].CloseFile();
    if (buffer[dim].OpenFileRead( dimLogs[dim]  )) return -1;
    //mprintf("\tOpened %s\n", logname.c_str());
    // Read the remlog header.
    int numexchg = ReadRemlogHeader(buffer[dim], log_type, buffer.size());
    if (numexchg == -1) return -1;
    mprintf("\t%s, type %s, should contain %i exchanges\n",
            dimLogs[dim].c_str(), LogDescription[log_type], numexchg);
    if (total_exchanges == -1)
      total_exchanges = numexchg;
    else if (numexchg != total_exchanges) {
      mprinterr("Error: Number of expected exchanges in dimension %i does not match\n"
                "Error: number of expected exchanges in first dimension.\n", dim + 1);
      return -1;
    }
    if (log_type != expectedType) {
      mprinterr("Error: Log '%s' is type %s, expected %s\n", dimLogs[dim].c_str(),
                LogDescription[log_type], LogDescription[expectedType]);
      return -1;
    }
  }
  return total_exchanges;
}

// -----------------------------------------------------------------------------
// DataIO_RemLog::ReadHelp()
void DataIO_RemLog::ReadHelp() {
  mprintf("\tnosearch            : Do not automatically search for MREMD dimension logs.\n"
          "\tdimfile <file>      : remd.dim file for processing MREMD logs.\n"
          "\tcrdidx <crd indices>: Use comma-separated list of indices as the initial\n"
          "\t                      coordinate indices.\n"
          "\tMultiple REM logs may be specified.\n");
}

// DataIO_RemLog::processReadArgs()
int DataIO_RemLog::processReadArgs(ArgList& argIn) {
  searchForLogs_ = !argIn.hasKey("nosearch");
  // Get dimfile arg
  dimfile_ = argIn.GetStringKey("dimfile");
  crdidx_ = argIn.GetStringKey("crdidx");
  // Unlike other data sets, remlog will find all file names now
  // in case of MREMD. Always add at least one entry to this array
  // for the first file.
  logFilenames_.push_back( std::string("") );
  std::string log_name = argIn.GetStringNext();
  while (!log_name.empty()) {
    FileName log(log_name);
    if (!File::Exists( log ))
      File::ErrorMsg( log.full() );
    else
      logFilenames_.push_back( log.Full() );
    log_name = argIn.GetStringNext();
  }
  return 0;
}

/// Get filename up to extension
//TODO May not need to be its own function. Make a general FileName function?
static inline std::string GetPrefix(FileName const& fname) {
  size_t found = fname.Full().rfind( fname.Ext() );
  return fname.Full().substr(0, found);
}

// DataIO_RemLog::ReadData()
int DataIO_RemLog::ReadData(FileName const& fnameIn, 
                            DataSetList& datasetlist, std::string const& dsname)
{
  bool processMREMD_ = false;
  if (!File::Exists( fnameIn )) {
    File::ErrorMsg( fnameIn.full() );
    return 1;
  }
  if (logFilenames_.empty()) // processReadArgs not called
    logFilenames_.push_back( fnameIn.Full() );
  else
    logFilenames_[0] = fnameIn.Full();
  DataSet_RemLog::GdimArray GroupDims; // TODO should only be in the DataSet
  ReplicaDimArray DimTypes;            // TODO should only be in the DataSet
  // Total number of replicas
  int n_mremd_replicas = 0;
  // Read replica dimension file if specified. Implies MREMD.
  if (!dimfile_.empty()) {
    // Sets up groups and dimensions
    n_mremd_replicas = ReadRemdDimFile( dimfile_, GroupDims, DimTypes );
    if (n_mremd_replicas < 1) {
      mprinterr("Error: Reading remd.dim file '%s'\n", dimfile_.c_str());
      return 1;
    }
    mprintf("\tExpecting %zu replica dimensions.\n", GroupDims.size());
  }
  // Split up crdidx arg
  ArgList idxArgs( crdidx_, "," );
  mprintf("\tReading from log files:");
  for (Sarray::const_iterator it = logFilenames_.begin(); it != logFilenames_.end(); ++it)
    mprintf(" %s", it->c_str());
  mprintf("\n");

  // Expected type for all logs
  LogType log_type;
  // Temperature map for dimensions (if needed) 
  std::vector<TmapType> TemperatureMap;
  // Coordinate indices for temperature dimensions (if needed)
  std::vector< std::vector<int> > TempCrdIdxs;
  // Log files for each dimension
  std::vector<BufferedLine> buffer( GroupDims.size() );
  // logFileGroups will be used to hold all dimension replica logs for all runs.
  //   When there are multiple dimensions, logFileGroups will look like:
  // {Run1Log.Dim1, Run1Log.Dim2, ... Run1Log.DimN}, {Run2Log.Dim1 ...}, ...
  //   otherwise:
  // {Run1Log}, {Run2Log}, ...
  typedef std::vector<Sarray> LogGroupType;
  LogGroupType logFileGroups;
  if (GroupDims.empty()) {
    // -------------------------------------------
    // Not M-REMD; T-REMD or H-REMD. Single dimension.
    processMREMD_ = false;
    for (Sarray::const_iterator logfile = logFilenames_.begin();
                                logfile != logFilenames_.end(); ++logfile)
      logFileGroups.push_back( Sarray( 1, *logfile ) );
    // Open the first log file specified to determine type.
    TemperatureMap.resize( 1 );
    TempCrdIdxs.resize( 1 );
    buffer.resize( 1 );
    if (buffer[0].OpenFileRead( logFilenames_.front() )) return 1;
    int numexchg = ReadRemlogHeader( buffer[0], log_type, 1 );
    if (numexchg < 1) return 1;
    ReplicaDimArray::RemDimType exch_type = ReplicaDimArray::UNKNOWN;
    if      (log_type == TREMD ) exch_type = ReplicaDimArray::TEMPERATURE;
    else if (log_type == HREMD ) exch_type = ReplicaDimArray::HAMILTONIAN;
    else if (log_type == RXSGLD) exch_type = ReplicaDimArray::RXSGLD;
    else if (log_type == PHREMD) exch_type = ReplicaDimArray::PH;
    else {
      mprinterr("Error: Invalid log type for single dimension: %s\n", LogDescription[log_type]);
      return 1;
    }
    DimTypes.AddRemdDimension( exch_type );
    // Set up dimension
    int group_size = 0;
    switch (exch_type) {
      case ReplicaDimArray::TEMPERATURE:
        TemperatureMap[0] = SetupTemperatureMap( buffer[0], TempCrdIdxs[0] );
        group_size = (int)TemperatureMap[0].size();
        mprintf("\t\t%i Temperature reps.\n", group_size);
        break;
      case ReplicaDimArray::PH:
        TemperatureMap[0] = Setup_pH_Map( buffer[0], TempCrdIdxs[0] );
        group_size = (int)TemperatureMap[0].size();
        mprintf("\t\t%i pH reps.\n", group_size);
        break;
      case ReplicaDimArray::HAMILTONIAN:
        group_size = CountHamiltonianReps( buffer[0] );
        mprintf("\t\t%i Hamiltonian reps.\n", group_size);
        break;
      case ReplicaDimArray::RXSGLD:
        group_size = CountHamiltonianReps( buffer[0] );
        mprintf("\t\t%i RXSGLD reps.\n", group_size);
        break;
      default: return 1; // sanity check
    }
    buffer[0].CloseFile();
    // In 1D number of replicas is just the group size
    n_mremd_replicas = group_size;
  } else {
    // -------------------------------------------
    // M-REMD
    processMREMD_ = true;
    // Ensure that each replica log has counterparts for each dimension
    // TODO: Read all headers and check dimensions in log.
    // Two cases; all log files were specified, or only lowest logs were specified.
    if ( searchForLogs_ ) { 
      FileName fname;
      for (Sarray::const_iterator logfile = logFilenames_.begin();
                                  logfile != logFilenames_.end(); ++logfile)
      {
        Sarray dimLogs;
        fname.SetFileName( *logfile );
        // Remove leading '.'
        std::string logExt = fname.Ext();
        if (logExt[0] == '.') logExt.erase(0,1);
        if ( !validInteger(logExt) ) {
          mprinterr("Error: MREMD log %s does not have valid numerical extension.\n", fname.full());
          return 1;
        }
        std::string Prefix = GetPrefix( fname );
        int numericalExt = convertToInteger( logExt );
        if (numericalExt != 1) {
          mprinterr("Error: Must specify MREMD log for dimension 1 (i.e. '%s.1')\n", 
                    Prefix.c_str());
          return 1;
        }
        dimLogs.push_back( *logfile );
        for (int idim = 2; idim <= (int)GroupDims.size(); idim++) {
          std::string logname = Prefix + "." + integerToString( idim );
          if ( !File::Exists(logname) ) {
            mprinterr("Error: MREMD log not found for dimension %i, '%s'\n",
                      idim, logname.c_str());
            return 1;
          }
          dimLogs.push_back( logname );
        }
        logFileGroups.push_back( dimLogs );
      }
    } else {
      // All logs specified. Assume they are given in order.
      Sarray dimLogs;
      Sarray::const_iterator logfile = logFilenames_.begin();
      while (logfile != logFilenames_.end()) {
        dimLogs.clear();
        for (unsigned int dim = 0; dim < GroupDims.size(); dim++) {
          if (logfile == logFilenames_.end()) {
            mprinterr("Error: Ran out of MREMD logs, run %zu, dimension %u\n",
                      logFileGroups.size() + 1, dim + 1);
            return 1;
          }
          dimLogs.push_back( *(logfile++) );
        }
        logFileGroups.push_back( dimLogs );
      }
    }
    // Set up temperature maps/coordinate index arrays for each dim/group.
    // Base this on the first set of MREMD replica logs.
    // Open remlogs for each dimension as buffered file.
    TemperatureMap.resize( GroupDims.size() );
    TempCrdIdxs.resize( GroupDims.size() );
    buffer.resize( GroupDims.size() );
    log_type = MREMD;
    int total_exchanges = OpenMremdDims(buffer, logFileGroups.front(), log_type);
    if (total_exchanges == -1) return 1;
    // Should now be positioned at the first exchange in each dimension.
    // Set up map/coordinate indices for each group and make sure they match
    // whats in the remd.dim file.
    for (int dim = 0; dim < (int)GroupDims.size(); dim++) {
      int group_size = 0;
      if ( DimTypes[dim] == ReplicaDimArray::TEMPERATURE ) {
        TemperatureMap[dim] = SetupTemperatureMap( buffer[dim], TempCrdIdxs[dim] );
        group_size = (int)TemperatureMap[dim].size();
        mprintf("\t\tDim %i: %i Temperature reps.\n", dim+1, group_size);
      } else if (DimTypes[dim] == ReplicaDimArray::PH) {
        TemperatureMap[dim] = Setup_pH_Map( buffer[dim], TempCrdIdxs[dim] );
        group_size = (int)TemperatureMap[dim].size();
        mprintf("\t\tDim %i: %i pH reps.\n", dim+1, group_size);
      } else if (DimTypes[dim] == ReplicaDimArray::HAMILTONIAN) {
        group_size = CountHamiltonianReps( buffer[dim] );
        mprintf("\t\tDim %i: %i Hamiltonian reps.\n", dim+1, group_size);
      } else if (DimTypes[dim] == ReplicaDimArray::RXSGLD) {
        group_size = CountHamiltonianReps( buffer[dim] );
        mprintf("\t\tDim %i: %i RXSGLD reps.\n", dim+1, group_size);
      } else {
        mprinterr("Error: Unrecognized dimension type.\n");
        return 1;
      }
      for (unsigned int grp = 0; grp < GroupDims[dim].size(); grp++) {
        if (group_size != (int)GroupDims[dim][grp].size()) {
          mprinterr("Error: Number of replicas in dimension %i (%i) does not match\n"
                    "Error: number of replicas in remd.dim file (%u)\n",
                    dim+1, group_size, GroupDims[dim][grp].size());
          return 1;
        }
      }
    } // END loop over replica dimensions
  }
  // SANITY CHECK
  if (n_mremd_replicas < 1) {
    mprinterr("Error: No replicas detected.\n");
    return 1;
  }

  // Coordinate indices for each replica. Start crdidx = repidx (from 1) for now.
  std::vector<int> CoordinateIndices( n_mremd_replicas );
  for (int repidx = 0; repidx < n_mremd_replicas; repidx++) {
    if (!idxArgs.empty()) {
      // User-specified starting coord indices
      CoordinateIndices[repidx] = idxArgs.getNextInteger(0);
      if (CoordinateIndices[repidx] < 0 || CoordinateIndices[repidx] > n_mremd_replicas )
      {
        mprinterr("Error: Given coordinate index out of range or not enough indices given.\n");
        return 1;
      }
    } else if (!processMREMD_ && !TempCrdIdxs.front().empty())
      // Use coordinate indices from 1D T-REMD log
      CoordinateIndices[repidx] = TempCrdIdxs.front()[repidx];
    else
      // Default: starting crdidx = repidx
      CoordinateIndices[repidx] = repidx + 1;
  }
  // Allocate replica log DataSet
  DataSet* ds = 0;
  if (!dsname.empty()) ds = datasetlist.CheckForSet( dsname );
  if (ds == 0) {
    // New set
    ds = datasetlist.AddSet( DataSet::REMLOG, dsname, "remlog" );
    if (ds == 0) return 1;
    ((DataSet_RemLog*)ds)->AllocateReplicas(n_mremd_replicas, GroupDims, DimTypes, 1, true, debug_);
  } else {
    if (ds->Type() != DataSet::REMLOG) {
      mprinterr("Error: Set '%s' is not replica log data.\n", ds->legend());
      return 1;
    }
    if ((int)ds->Size() != n_mremd_replicas) {
      mprinterr("Error: Replica log data '%s' is set up for %zu replicas,"
                " current # replicas is %i\n", ds->legend(), ds->Size(),
                n_mremd_replicas);
      return 1;
    }
    mprintf("\tReading final coordinate indices from last frame of existing set.\n");
    for (int repidx = 0; repidx < n_mremd_replicas; repidx++)
      CoordinateIndices[repidx] = ((DataSet_RemLog*)ds)->LastRepFrame(repidx).CoordsIdx();
  }
//  if (!idxArgs.empty()) {
    mprintf("\tInitial coordinate indices:");
    for (std::vector<int>::const_iterator c = CoordinateIndices.begin();
                                          c != CoordinateIndices.end(); ++c)
      mprintf(" %i", *c);
    mprintf("\n");
//  }
  // Loop over all remlogs
  DataSet_RemLog& ensemble = static_cast<DataSet_RemLog&>( *ds );
  for (LogGroupType::const_iterator it = logFileGroups.begin(); it != logFileGroups.end(); ++it)
  { 
    // Open the current remlog, advance to first exchange
    int numexchg = OpenMremdDims(buffer, *it, log_type);
    if (numexchg == -1) return 1;
    mprintf("\t%s should contain %i exchanges\n", it->front().c_str(), numexchg);
    // Should now be positioned at 'exchange 1'.
    // Loop over all exchanges.
    ProgressBar progress( numexchg );
    bool fileEOF = false;
    const char* ptr = 0;
    unsigned int current_dim = 0;
    int grp; // Will be set to group number for MREMD or group index otherwise
    for (int exchg = 0; exchg < numexchg; exchg++) {
      progress.Update( exchg );
      // Loop over all groups in the current dimension
      for (unsigned int gidx = 0; gidx < ensemble.GroupDims()[current_dim].size(); gidx++)
      {
        if (processMREMD_) {
          if (sscanf(buffer[current_dim].CurrentLine(), "%*s%*s%*i%*s%*s%i", &grp)!=1) {
            mprinterr("Error: Could not get MREMD group number.\n");
            return 1;
          }
          // REMD group indices start from 1
          grp--;
        } else
          grp = gidx;
        // Loop over all replicas in the current group
        //mprintf("--------------------------------------------------------------------------\n");
        for (unsigned int replica = 0;
                          replica < ensemble.GroupDims()[current_dim][grp].size(); replica++)
        {
          // Read remlog line.
          ptr = buffer[current_dim].Line();
          if (ptr == 0) {
            mprinterr("Error: reading remlog; unexpected EOF. Dim=%u, Exchg=%i, grp=%u, rep=%u\n",
                      current_dim+1, exchg+1, grp+1, replica+1);
            fileEOF = true;
            // If this is not the first replica remove all partial replicas
            if (replica > 0) ensemble.TrimLastExchange();
            break;
          }
          // ----- T-REMD ----------------------------
          /* Format:
           * '(i2,6f10.2,i8)'
          # Rep#, Velocity Scaling, T, Eptot, Temp0, NewTemp0, Success rate (i,i+1), ResStruct#
            1     -1.00      0.00   -433.24    300.00    300.00      0.00      -1
           * Order during REMD is exchange -> MD, so NewTemp0 is the temp. that gets
           * simulated. TODO: Is that valid?
           */
          if (DimTypes[current_dim] == ReplicaDimArray::TEMPERATURE) {
            int tremd_crdidx, current_crdidx; // TODO: Remove tremd_crdidx
            double tremd_scaling, tremd_pe, tremd_temp0, tremd_tempP;
            if ( sscanf(ptr, "%2i%10lf%*10f%10lf%10lf%10lf", &tremd_crdidx, &tremd_scaling,
                        &tremd_pe, &tremd_temp0, &tremd_tempP) != 5 )
            {
              mprinterr("Error reading TREMD line from rem log. Dim=%u, Exchg=%i, grp=%u, rep=%u\n",
                        current_dim+1, exchg+1, grp+1, replica+1);
              mprinterr("Error: Line: %s", ptr);
              return 1;
            }
            // Figure out my position within the group.
            TmapType::const_iterator tmap = 
              TemperatureMap[current_dim].find( tremd_temp0 );
            if (tmap == TemperatureMap[current_dim].end()) {
              mprinterr("Error: replica temperature %.2f not found in temperature map.\n", 
                        tremd_temp0);
              return 1;
            }
            // What is my actual position? Currently mapped rep nums start from 1
            int tremd_repidx = ensemble.GroupDims()[current_dim][grp][tmap->second - 1].Me();
            // Who is my partner? ONLY VALID IF EXCHANGE OCCURS
            tmap = TemperatureMap[current_dim].find( tremd_tempP ); // TODO: Make function
            if (tmap == TemperatureMap[current_dim].end()) {
              mprinterr("Error: partner temperature %.2f not found in temperature map.\n",
                        tremd_tempP);
              return 1;
            }
            int tremd_partneridx = ensemble.GroupDims()[current_dim][grp][tmap->second - 1].Me();
            // Exchange success if velocity scaling is > 0.0
            bool tremd_success = (tremd_scaling > 0.0);
            // If an exchange occured, coordsIdx will be that of partner replica.
            if (tremd_success)
              current_crdidx = CoordinateIndices[tremd_partneridx-1];
            else
              current_crdidx = CoordinateIndices[tremd_repidx-1];
            //mprintf("DEBUG: Exchg %8i Tdim# %2u T=%6.2f group=%2u"
            //        " repidx=%3i partneridx=%3i oldcrdidx=%i newcrdidx=%i\n",
            //        exchg+1, current_dim+1, tremd_temp0, grp+1,
            //        tremd_repidx, tremd_partneridx,
            //        CoordinateIndices[tremd_repidx-1], current_crdidx);
            // Create replica frame for TREMD
            ensemble.AddRepFrame( tremd_repidx-1,
                                  DataSet_RemLog:: 
                                  ReplicaFrame(tremd_repidx, tremd_partneridx,
                                               current_crdidx, current_dim,
                                               tremd_success,
                                               tremd_temp0, tremd_pe, 0.0) );
          // ----- pH-REMD ----------------------------
          /* Format:
           * '(i6,x,i7,x,2f7.3,x,f8.4)'
           # Rep#, N_prot, old_pH, new_pH, Success rate (i,i+1)
           # exchange        1
                1       1   2.000  3.500   0.0000
           * Order during REMD is exchange -> MD, so old_pH is the pH that gets
           * simulated. TODO: Is that valid?
           */
          } else if (DimTypes[current_dim] == ReplicaDimArray::PH) {
            int ph_crdidx, current_crdidx; // TODO: Remove ph_crdidx
            double old_pH, new_pH;
            if ( sscanf(ptr, "%6i%*8i%*c%7lf%7lf", &ph_crdidx, &old_pH, &new_pH) != 3)
            {
              mprinterr("Error reading PH line from rem log. Dim=%u, Exchg=%i, grp=%u, rep=%u\n",
                        current_dim+1, exchg+1, grp+1, replica+1);
              mprinterr("Error: Line: %s", ptr);
              return 1;
            }
            // Figure out my position within the group.
            TmapType::const_iterator pmap = 
              TemperatureMap[current_dim].find( old_pH );
            if (pmap == TemperatureMap[current_dim].end()) {
              mprinterr("Error: replica pH %.2f not found in pH map.\n", old_pH);
              return 1;
            }
            // What is my actual position? Currently mapped rep nums start from 1
            int phremd_repidx = ensemble.GroupDims()[current_dim][grp][pmap->second - 1].Me();
            // Who is my partner? ONLY VALID IF EXCHANGE OCCURS
            pmap = TemperatureMap[current_dim].find( new_pH ); // TODO: Make function
            if (pmap == TemperatureMap[current_dim].end()) {
              mprinterr("Error: partner pH %.2f not found in pH map.\n", new_pH);
              return 1;
            }
            int phremd_partneridx = ensemble.GroupDims()[current_dim][grp][pmap->second - 1].Me();
            // Exchange success if old_pH != new_pH
            double delta_pH = new_pH - old_pH;
            bool phremd_success = (delta_pH > 0.0 || delta_pH < 0.0);
            // If an exchange occured, coordsIdx will be that of partner replica.
            if (phremd_success)
              current_crdidx = CoordinateIndices[phremd_partneridx-1];
            else
              current_crdidx = CoordinateIndices[phremd_repidx-1];
            //mprintf("DEBUG: Exchg %3i pHdim# %1u old_pH=%6.2f group=%2u"
            //        " repidx=%3i partneridx=%3i oldcrdidx=%i newcrdidx=%i Exch=%i\n",
            //        exchg+1, current_dim+1, old_pH, grp+1,
            //        phremd_repidx, phremd_partneridx,
            //        CoordinateIndices[phremd_repidx-1], current_crdidx, (int)phremd_success);
            // Create replica frame for PHREMD
            ensemble.AddRepFrame( phremd_repidx-1,
                                  DataSet_RemLog:: 
                                  ReplicaFrame(phremd_repidx, phremd_partneridx,
                                               current_crdidx, current_dim,
                                               phremd_success,
                                               old_pH, 0.0, 0.0) );
          // ----- H-REMD ----------------------------
          /* Format:
           * '(2i6,5f10.2,4x,a,2x,f10.2)'
       # Rep#, Neibr#, Temp0, PotE(x_1), PotE(x_2), left_fe, right_fe, Success, Success rate (i,i+1)
            1     8    300.00 -25011.03 -24959.58    -27.48      0.00    F        0.00
           */
          } else if (DimTypes[current_dim] == ReplicaDimArray::HAMILTONIAN) {
            int hremd_grp_repidx, hremd_grp_partneridx, current_crdidx;
            double hremd_temp0, hremd_pe_x1, hremd_pe_x2;
            bool hremd_success;
            if ( sscanf(ptr, "%6i%6i%10lf%10lf%10lf", &hremd_grp_repidx, &hremd_grp_partneridx,
                        &hremd_temp0, &hremd_pe_x1, &hremd_pe_x2) != 5 )
            {
              mprinterr("Error reading HREMD line from rem log. Dim=%u, Exchg=%i, grp=%u, rep=%u\n",
                        current_dim+1, exchg+1, grp+1, replica+1);
              mprinterr("Error: Line: %s", ptr);
              return 1;
            }
            // What is my actual position and who is my actual partner?
            int hremd_repidx = ensemble.GroupDims()[current_dim][grp][hremd_grp_repidx-1].Me();
            int hremd_partneridx = ensemble.GroupDims()[current_dim][grp][hremd_grp_partneridx-1].Me();
            //mprintf("DEBUG: Exchg %i Hdim# %u group=%u group_repidx=%i repidx=%i\n",
            //        exchg+1, current_dim+1, grp+1, hremd_grp_repidx, hremd_repidx);
            // Determine if an exchange occurred
            switch ( ptr[66] ) {
              case 'T': hremd_success = true; break;
              case 'F': hremd_success = false; break;
              default: // Should only get here with malformed HREMD log file.
                mprinterr("Error: expected only 'T' or 'F' at character 67, got %c\n", ptr[66]);
                return 1;
            }
            // If an exchange occured, coordsIdx will be that of partner replica.
            if (hremd_success)
              current_crdidx = CoordinateIndices[hremd_partneridx-1];
            else
              current_crdidx = CoordinateIndices[hremd_repidx-1];
            // Create replica frame for HREMD
            ensemble.AddRepFrame( hremd_repidx-1,
                                  DataSet_RemLog::
                                  ReplicaFrame(hremd_repidx, hremd_partneridx,
                                               current_crdidx, current_dim,
                                               hremd_success,
                                               hremd_temp0, hremd_pe_x1, hremd_pe_x2) );
          // ----- RXSGLD ----------------------------
          /* Format:
           * (i4,i4,2f8.4,2f8.2,e14.6,f8.4) 
           # Rep Stagid Vscale  SGscale Temp Templf Eptot Acceptance(i,i+1)
              1   1  1.0000  1.0000    0.00   21.21 -0.102302E+02  0.0000
           */
          } else if (DimTypes[current_dim] == ReplicaDimArray::RXSGLD) {
            // Consider accept if sgscale is not -1.0.
            int sgld_repidx, sgld_crdidx;
            double sgscale;
            if ( sscanf(ptr, "%4i%4i%*8f%8lf", &sgld_crdidx, &sgld_repidx, &sgscale) != 3 ) {
              mprinterr("Error reading RXSGLD line from rem log. "
                        "Dim=%u, Exchg=%i, grp=%u, rep=%u\n",
                        current_dim+1, exchg+1, grp+1, replica+1);
              mprinterr("Error: Line: %s", ptr);
              return 1;
            }
            bool sgld_success = (sgscale > -1.0);
            // The partner index is not stored in RXSGLD logs.
            bool sgld_up;
            // First exchange for even replicas is up, then down.
            if ((exchg % 2)==0) {
              if ((sgld_repidx % 2)==0) // Exchange up
                sgld_up = true;
              else
                sgld_up = false;
            } else {
              if ((sgld_repidx % 2)==0) // Exchange down
                sgld_up = false;
              else
                sgld_up = true;
            }
            int sgld_partneridx;
            if (sgld_up)
              sgld_partneridx = ensemble.GroupDims()[current_dim][grp][sgld_repidx-1].R_partner();
            else
              sgld_partneridx = ensemble.GroupDims()[current_dim][grp][sgld_repidx-1].L_partner();
            // If an exchange occured, coordsIdx will be that of partner replica.
            int current_crdidx;
            if (sgld_success)
              current_crdidx = CoordinateIndices[sgld_partneridx-1];
            else
              current_crdidx = CoordinateIndices[sgld_repidx-1];
            // Create replica frame for SGLD
            ensemble.AddRepFrame( sgld_repidx-1,
                                  DataSet_RemLog::
                                  ReplicaFrame(sgld_repidx, sgld_partneridx,
                                               current_crdidx, current_dim,
                                               sgld_success,
                                               0.0, 0.0, 0.0) );
          // -----------------------------------------
          } else {
            mprinterr("Error: remlog; unknown type.\n");
            return 1;
          }
          // -----------------------------------------
        } // END loop over replicas in group
        if ( fileEOF ) break; // Error occurred reading replicas, skip rest of groups.
        // Read next group exchange line.
        ptr = buffer[current_dim].Line();
      } // END loop over groups in dimension
      if ( fileEOF ) break; // Error occurred reading replicas, skip rest of exchanges.
      // Update coordinate indices.
      //mprintf("DEBUG: exchange= %i: Updating coordinates\n", exchg + 1);
      for (int repidx = 0; repidx < n_mremd_replicas; repidx++) {
        //mprintf("DEBUG:\tReplica %i crdidx %i =>", repidx+1, CoordinateIndices[repidx]);
        CoordinateIndices[repidx] = ensemble.LastRepFrame(repidx).CoordsIdx();
        //mprintf(" %i\n", CoordinateIndices[repidx]); // DEBUG
      }
      // Currently each exchange the dimension alternates
      ++current_dim;
      if (current_dim == ensemble.GroupDims().size()) current_dim = 0;
    } // END loop over exchanges in remlog
  } // END loop over remlogs
  if (!ensemble.ValidEnsemble()) {
    mprinterr("Error: Ensemble is not valid.\n");
    return 1;
  }
  if (debug_ > 1)
    ensemble.PrintReplicaStats();
  // ---------------------------------------------

  return 0;
}
