#include <cstdio> // sscanf
#include "DataIO_RemLog.h"
#include "CpptrajStdio.h" 
#include "ProgressBar.h"
#include "DataSet_RemLog.h"

// CONSTRUCTOR
DataIO_RemLog::DataIO_RemLog() : debug_(0) {}

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

// DataIO_RemLog::ReadData()
int DataIO_RemLog::ReadData(std::string const& fname, ArgList& argIn,
                            DataSetList& datasetlist, std::string const& dsname)
{
  ExchgType firstlog_type = UNKNOWN;
  std::vector<std::string> logFilenames;
  logFilenames.push_back( fname );
  // Check if more than one log name was specified.
  ArgList lognames(argIn.GetStringKey("logfiles"), ",");
  if (!lognames.empty()) {
    for (int i = 0; i < lognames.Nargs(); i++)
      logFilenames.push_back( lognames[i] );
  }
  mprintf("\tReading from log files:");
  for (std::vector<std::string>::const_iterator it = logFilenames.begin();
                                          it != logFilenames.end(); ++it)
    mprintf(" %s", (*it).c_str());
  mprintf("\n");
  // Open first remlog as buffered file
  BufferedLine buffer;
  if (buffer.OpenFileRead( fname )) return 1;
  // Read the first line. Should be '# Replica Exchange log file'
  if (ReadRemlogHeader(buffer, firstlog_type) == -1) return 1;
  // Should currently be positioned at the first exchange. Need to read this
  // to determine how many replicas there are (and temperature for T-REMD).
  DataSet_RemLog::TmapType TemperatureMap;
  double t0;
  int n_replicas = 0;
  const char* ptr = buffer.Line();
  while (ptr != 0 && ptr[0] != '#') {
    if (firstlog_type == TREMD) {
      // For temperature remlog create temperature map. Map temperatures to 
      // index + 1 since indices in the remlog start from 1.
      //mprintf("DEBUG: Temp0= %s", ptr+32);
      if ( sscanf(ptr+32, "%10lf", &t0) != 1) {
        mprinterr("Error: could not read temperature from T-REMD log.\n"
                  "Error: Line: %s", ptr);
        return 1;
      }
      std::pair<DataSet_RemLog::TmapType::iterator,bool> ret =
        TemperatureMap.insert(std::pair<double,int>( t0, n_replicas+1 ));
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
  DataSet* ds = datasetlist.AddSet( DataSet::REMLOG, dsname, "remlog" );
  if (ds == 0) return 1;
  DataSet_RemLog& ensemble = static_cast<DataSet_RemLog&>( *ds );
  ensemble.AllocateReplicas(n_replicas);
  std::vector<DataSet_RemLog::ReplicaFrame> replicaFrames;
  std::vector<int> coordinateIndices;
  if (firstlog_type == HREMD) {
    replicaFrames.resize( n_replicas );
    // All coord indices start equal to replica indices.
    // Indices start from 1 in remlogs (H-REMD only).
    coordinateIndices.resize( n_replicas );
    for (int replica = 0; replica < n_replicas; replica++)
      coordinateIndices[replica] = replica+1;
  } else if (firstlog_type == TREMD)
    replicaFrames.resize(1);
  // Close first remlog 
  buffer.CloseFile();

  for (std::vector<std::string>::const_iterator it = logFilenames.begin();
                                                it != logFilenames.end(); ++it)
  {
    // Open the current remlog, advance to first exchange
    buffer.OpenFileRead( *it );
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
          if (replicaFrames[replica].SetHremdFrame( ptr, coordinateIndices[replica] )) {
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
        // Determine whether exchanges occurred. Update coordinate indices accordingly.
        //mprintf("DEBUG: exchange= %i:\n", exchg + 1);
        for (int replica = 0; replica < n_replicas; replica++) {
          //mprintf("DEBUG:\tReplica %i crdidx %i =>", replica+1, coordinateIndices[replica]);
          if (replicaFrames[replica].Success()) {
            int partner = replicaFrames[replica].PartnerIdx() - 1;
            coordinateIndices[replica] = replicaFrames[partner].CoordsIdx();
          }
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
    for (int j = 0; j < ensemble.Size(); j++) {
      DataSet_RemLog::ReplicaFrame const& frm = ensemble.RepFrame(i, j); 
      mprintf("%10u %6i %6i %6i %12.4f %12.4f %12.4f %1i\n", i + 1,
              frm.ReplicaIdx(), frm.PartnerIdx(), frm.CoordsIdx(), frm.Temp0(), 
              frm.PE_X1(), frm.PE_X2(), (int)frm.Success());
    }
  } 
*/
  return 0;
}
