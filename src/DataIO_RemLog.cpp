#include <cstdio> // sscanf
#include <set> // for TREMD temperature sorting
#include "DataIO_RemLog.h"
#include "CpptrajStdio.h" 
#include "ProgressBar.h"
#include "DataSet_RemLog.h"
#include "StringRoutines.h" // fileExists

// CONSTRUCTOR
DataIO_RemLog::DataIO_RemLog() : debug_(0) {
  SetValid( DataSet::REMLOG );
}

const char* ExchgDescription[] = {
"Unknown", "Temperature", "Hamiltonian", "MultipleDim"
};

// DataIO_RemLog::ReadRemlogHeader()
int DataIO_RemLog::ReadRemlogHeader(BufferedLine& buffer, ExchgType& type) {
  int numexchg = -1;
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
    if (columns.hasKey("exchange")) break;
    if (debug_ > 0) mprintf("\t%s", line.c_str());
    if (columns.hasKey("numexchg")) {
      numexchg = columns.getNextInteger(-1);
    }
    if (columns.hasKey("Rep#,")) {
      if (columns[2] == "Neibr#,") type = HREMD;
      else if (columns[2] == "Velocity") type = TREMD;
    }
  }
  if (numexchg < 1) {
    mprinterr("Error: Invalid number of exchanges (%i) in rem log.\n");
    return -1;
  }
  return numexchg;
}

static inline bool IsNullPtr( const char* ptr ) { 
  if (ptr == 0) {
    mprinterr("Error: Could not read file.\n");
    return true;
  }
  return false;
}

static const char* dimTypeString[] = {0, "TREMD", "HREMD", "MREMD"};

// DataIO_RemLog::ReadRemdDimFile()
int DataIO_RemLog::ReadRemdDimFile(std::string const& rd_name) {
  BufferedLine rd_file;
  if (rd_file.OpenFileRead( rd_name )) {
    mprinterr("Error: Could not read remd dim file '%s'\n", rd_name.c_str());
    return 1;
  }
  const char* separators = " =,()";
  // Read dimension title
  const char* ptr = rd_file.Line();
  if (IsNullPtr( ptr )) return 1;
  // ptr Should end with a newline
  mprintf("\tReplica dimension file '%s' title: %s", rd_name.c_str(), ptr);
  // Read each &multirem section
  GroupDims_.clear();
  ArgList rd_arg;
  while (ptr != 0) {
    rd_arg.SetList( std::string(ptr), separators );
    if ( rd_arg[0] == "&multirem" ) {
      GroupDimType Groups;
      std::string desc;
      ExchgType exch_type = UNKNOWN;
      while (ptr != 0) {
        rd_arg.SetList( std::string(ptr), separators );
        if (rd_arg.CommandIs("&end") || rd_arg.CommandIs("/")) break;
        rd_arg.MarkArg(0);
        if ( rd_arg.CommandIs("exch_type") ) {
          if ( rd_arg.hasKey("TEMP") || rd_arg.hasKey("TEMPERATURE") )
            exch_type = TREMD;
          else if ( rd_arg.hasKey("HAMILTONIAN") || rd_arg.hasKey("HREMD") )
            exch_type = HREMD;
          else {
            mprinterr("Error: Unrecognized exch_type: %s\n", rd_arg.ArgLine());
            return 1;
          }
        } else if ( rd_arg.CommandIs("group") ) {
          int group_num = rd_arg.getNextInteger(-1);
          if (group_num < 1) {
            mprinterr("Error: Invalid group number: %i\n", group_num);
            return 1;
          }
          mprintf("\t\tGroup %i\n", group_num);
          std::vector<int> indices;
          int group_index = rd_arg.getNextInteger(-1);
          while (group_index != -1) {
            indices.push_back( group_index );
            group_index = rd_arg.getNextInteger(-1);
          }
          // Set up partner array for this group
          GroupArray group;
          for (int i = 0; i < (int)indices.size(); i++) {
            int l_idx = i - 1;
            if (l_idx < 0) l_idx = (int)indices.size() - 1;
            int r_idx = i + 1;
            if (r_idx == (int)indices.size()) r_idx = 0;
            group.push_back( GroupReplica(indices[l_idx], indices[i], indices[r_idx]) );
            mprintf("\t\t\t%i - %i - %i\n", group.back().L_partner(),
                    group.back().Me(), group.back().R_partner());
          }
          Groups.push_back( group ); 
        } else if ( rd_arg.CommandIs("desc") ) {
          desc = rd_arg.GetStringNext(); 
        }
        ptr = rd_file.Line();
      }
      mprintf("\tDimension %zu: type '%s' description '%s' groups=%zu\n", GroupDims_.size() + 1,
              dimTypeString[exch_type], desc.c_str(), Groups.size());
      GroupDims_.push_back( Groups ); 
    }
    ptr = rd_file.Line();
  }
  if (GroupDims_.empty()) {
    mprinterr("Error: No replica dimensions found.\n");
    return 1;
  }

  return 0;
}

void DataIO_RemLog::ReadHelp() {
  mprintf("\tcrdidx <crd indices>: Use comma-separated list of indices as the initial\n"
          "\t                      coordinate indices (H-REMD only).\n"
          "\tMultiple REM logs may be specified.\n");
}

static inline std::string GetPrefix(FileName const& fname) {
  size_t found = fname.Full().rfind( fname.Ext() );
  return fname.Full().substr(0, found);
}

// DataIO_RemLog::MremdRead()
int DataIO_RemLog::MremdRead(std::vector<std::string> const& logFilenames,
                             DataSetList& datasetlist, std::string const& dsname)
{
  // Ensure that each replica log has counterparts for each dimension
  FileName fname;
  for (std::vector<std::string>::const_iterator logfile = logFilenames.begin();
                                                logfile != logFilenames.end();
                                              ++logfile)
  {
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
      mprinterr("Error: Must specify MREMD log for dimension 1 (i.e. '%s.1')\n", Prefix.c_str());
      return 1;
    }
    for (int i = 2; i < (int)GroupDims_.size(); i++) {
      std::string logname = Prefix + "." + integerToString( i );
      if ( !fileExists(logname) ) {
        mprinterr("Error: MREMD log not found for dimension %i, '%s'\n",
                  i, logname.c_str());
        return 1;
      }
    }
  }

  return 0;
}
  
// DataIO_RemLog::ReadData()
int DataIO_RemLog::ReadData(std::string const& fname, ArgList& argIn,
                            DataSetList& datasetlist, std::string const& dsname)
{
  ExchgType firstlog_type = UNKNOWN;
  std::vector<std::string> logFilenames;
  if (!fileExists( fname )) {
    mprinterr("Error: File '%s' does not exist.\n", fname.c_str());
    return 1;
  }
  logFilenames.push_back( fname );
  // Get dimfile arg
  std::string dimfile = argIn.GetStringKey("dimfile");
  if (!dimfile.empty()) {
    if (ReadRemdDimFile( dimfile )) {
      mprinterr("Error: Reading remd.dim file '%s'\n", dimfile.c_str());
      return 1;
    }
    mprintf("\tExpecting %zu replica dimensions.\n", GroupDims_.size());
  }
  // Get crdidx arg
  ArgList idxArgs( argIn.GetStringKey("crdidx"), "," );
  // Check if more than one log name was specified.
  std::string log_name = argIn.GetStringNext();
  while (!log_name.empty()) {
    if (!fileExists( log_name ))
      mprintf("Warning: '%s' does not exist.\n", log_name.c_str());
    else
      logFilenames.push_back( log_name );
    log_name = argIn.GetStringNext();
  }
  mprintf("\tReading from log files:");
  for (std::vector<std::string>::const_iterator it = logFilenames.begin();
                                          it != logFilenames.end(); ++it)
    mprintf(" %s", (*it).c_str());
  mprintf("\n");
  // Multidim replica read requires reading from multiple files
  if (!GroupDims_.empty()) return MremdRead( logFilenames, datasetlist, dsname );
  // Open first remlog as buffered file
  BufferedLine buffer;
  if (buffer.OpenFileRead( fname )) return 1;
  // Read the first line. Should be '# Replica Exchange log file'
  if (ReadRemlogHeader(buffer, firstlog_type) == -1) return 1;
  // Should currently be positioned at the first exchange. Need to read this
  // to determine how many replicas there are (and temperature for T-REMD).
  DataSet_RemLog::TmapType TemperatureMap;
  std::set<double> tList;
  double t0;
  int n_replicas = 0;
  const char* ptr = buffer.Line();
  while (ptr != 0 && ptr[0] != '#') {
    if (firstlog_type == TREMD) {
      // For temperature remlog create temperature map.
      //mprintf("DEBUG: Temp0= %s", ptr+32);
      if ( sscanf(ptr+32, "%10lf", &t0) != 1) {
        mprinterr("Error: could not read temperature from T-REMD log.\n"
                  "Error: Line: %s", ptr);
        return 1;
      }
      std::pair<std::set<double>::iterator,bool> ret = tList.insert( t0 );
      if (!ret.second) {
        mprinterr("Error: duplicate temperature %.2f detected in T-REMD remlog\n", t0);
        return 1;
      }
    }
    ptr = buffer.Line();
    ++n_replicas;
  }
  mprintf("\t%i %s replicas.\n", n_replicas, ExchgDescription[firstlog_type]);
  if (n_replicas < 1) {
    mprinterr("Error: Detected less than 1 replica in remlog.\n");
    return 1;
  }
  if (firstlog_type == TREMD) {
    // Temperatures are already sorted lowest to highest in set. Map 
    // temperatures to index + 1 since indices in the remlog start from 1.
    int repnum = 1;
    for (std::set<double>::const_iterator temp0 = tList.begin(); temp0 != tList.end(); ++temp0)
      TemperatureMap.insert(std::pair<double,int>(*temp0, repnum++));
    for (DataSet_RemLog::TmapType::const_iterator tmap = TemperatureMap.begin();
                                                  tmap != TemperatureMap.end(); ++tmap)
      mprintf("\t\t%i => %f\n", (*tmap).second, (*tmap).first);
  }
  DataSet* ds = datasetlist.AddSet( DataSet::REMLOG, dsname, "remlog" );
  if (ds == 0) return 1;
  DataSet_RemLog& ensemble = static_cast<DataSet_RemLog&>( *ds );
  ensemble.AllocateReplicas(n_replicas);
  std::vector<DataSet_RemLog::ReplicaFrame> replicaFrames;
  std::vector<int> coordinateIndices;
  if (firstlog_type == HREMD) {
    replicaFrames.resize( n_replicas );
    if (!idxArgs.empty() && idxArgs.Nargs() != n_replicas) {
      mprinterr("Error: crdidx: Ensemble size is %i but only %i indices given!\n",
                n_replicas, idxArgs.Nargs());
      return 1;
    }
    // All coord indices start equal to replica indices.
    // Indices start from 1 in remlogs (H-REMD only).
    coordinateIndices.resize( n_replicas );
    mprintf("\tInitial H-REMD coordinate indices:");
    for (int replica = 0; replica < n_replicas; replica++) {
      if (idxArgs.empty())
        coordinateIndices[replica] = replica+1;
      else // TODO: Check replica index range
        coordinateIndices[replica] = idxArgs.getNextInteger(0);
      mprintf(" %i", coordinateIndices[replica]);
    }
    mprintf("\n");
  } else if (firstlog_type == TREMD)
    replicaFrames.resize(1);
  // Close first remlog 
  buffer.CloseFile();

  for (std::vector<std::string>::const_iterator it = logFilenames.begin();
                                                it != logFilenames.end(); ++it)
  {
    // Open the current remlog, advance to first exchange
    if (buffer.OpenFileRead( *it )) return 1;
    //ptr = buffer.Line();
    //while (ptr[0] == '#' && ptr[2] != 'e' && ptr[3] != 'x') ptr = buffer.Line();
    ExchgType thislog_type = UNKNOWN;
    int numexchg = ReadRemlogHeader(buffer, thislog_type);
    if (thislog_type != firstlog_type) {
      mprinterr("Error: rem log %s type %s does not match first rem log.\n",
                (*it).c_str(), ExchgDescription[thislog_type]);
      return 1;
    }
    mprintf("\t%s should contain %i exchanges\n", (*it).c_str(), numexchg);
    // Should now be positioned at 'exchange 1'.
    // Loop over all exchanges.
    ProgressBar progress( numexchg );
    bool fileEOF = false;
    for (int exchg = 0; exchg < numexchg; exchg++) {
      progress.Update( exchg );
      for (int replica = 0; replica < n_replicas; replica++) {
        // Read remlog line.
        ptr = buffer.Line();
        if (ptr == 0) {
          mprinterr("Error: reading remlog; unexpected end of file. Exchange=%i, replica=%i\n",
                    exchg+1, replica+1);
          fileEOF = true;
          // If this is not the first replica remove all partial replicas
          if (replica > 0) ensemble.TrimLastExchange();
          break;
        }
        // ----- T-REMD ----------------------------
        if (thislog_type == TREMD) {
          if (replicaFrames[0].SetTremdFrame( ptr, TemperatureMap )) {
          mprinterr("Error reading TREMD line from rem log. Exchange=%i, replica=%i\n",
                      exchg+1, replica+1);
            return 1;
          }
          // Add replica frame to appropriate ensemble
          ensemble.AddRepFrame( replicaFrames[0].ReplicaIdx()-1, replicaFrames[0] );
        // ----- H-REMD ----------------------------
        } else if (thislog_type == HREMD) {
          if (replicaFrames[replica].SetHremdFrame( ptr, coordinateIndices )) {
            mprinterr("Error reading HREMD line from rem log. Exchange=%i, replica=%i\n",
                      exchg+1, replica+1);
            return 1;
          }
          // Add replica frame to appropriate ensemble
          ensemble.AddRepFrame( replica, replicaFrames[replica] );
        // -----------------------------------------
        } else {
          mprinterr("Error: remlog; unknown type.\n");
        }
      }
      if ( fileEOF ) break; // Error occurred reading replicas, skip rest of exchanges.
      if (thislog_type == HREMD) {
        // Update coordinate indices.
        //mprintf("DEBUG: exchange= %i:\n", exchg + 1);
        for (int replica = 0; replica < n_replicas; replica++) {
          //mprintf("DEBUG:\tReplica %i crdidx %i =>", replica+1, coordinateIndices[replica]);
          coordinateIndices[replica] = replicaFrames[replica].CoordsIdx();
          //mprintf(" %i\n", coordinateIndices[replica]); // DEBUG
        }
      }
      // Read 'exchange N' line.
      ptr = buffer.Line();
    } // END loop over exchanges
    buffer.CloseFile();
  } // END loop over remlog files
  if (!ensemble.ValidEnsemble()) {
    mprinterr("Error: Ensemble is not valid.\n");
    return 1;
  }
  // DEBUG - Print out replica 1 stats
/*
  mprintf("Replica Stats:\n"
          "%-10s %6s %6s %6s %12s %12s %12s S\n", "#Exchange", "RepIdx", "PrtIdx", "CrdIdx",
          "Temp0", "PE_X1", "PE_X2");
  for (int i = 0; i < ensemble.NumExchange(); i++) {
    for (int j = 0; j < (int)ensemble.Size(); j++) {
      DataSet_RemLog::ReplicaFrame const& frm = ensemble.RepFrame(i, j); 
      mprintf("%10u %6i %6i %6i %12.4f %12.4f %12.4f %1i\n", i + 1,
              frm.ReplicaIdx(), frm.PartnerIdx(), frm.CoordsIdx(), frm.Temp0(), 
              frm.PE_X1(), frm.PE_X2(), (int)frm.Success());
    }
  } 
*/
  return 0;
}
