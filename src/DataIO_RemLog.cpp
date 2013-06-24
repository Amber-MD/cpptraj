#include <cstdio> // sscanf
#include "DataIO_RemLog.h"
#include "CpptrajStdio.h" 
#include "BufferedLine.h"
#include "ProgressBar.h"
#include "DataSet_RemLog.h"

// CONSTRUCTOR
DataIO_RemLog::DataIO_RemLog() {}

const char* ExchgDescription[] = {
"Unknown", "Temperature", "Hamiltonian", "MultipleDim"
};

// DataIO_RemLog::ReadData()
int DataIO_RemLog::ReadData(std::string const& fname, DataSetList& datasetlist) {
  enum ExchgType { UNKNOWN = 0, TREMD, HREMD, MREMD };
  ExchgType type = UNKNOWN;
  int numexchg = 0;
  int LineNum = 0;
  //const char* SEPARATORS = " ";
  // Buffer file
  BufferedLine buffer;
  if (buffer.OpenFileRead( fname )) return 1;

  // Read the first line. Should be '# Replica Exchange log file'
  std::string line = buffer.GetLine();
  ++LineNum;
  if (line.compare(0, 27, "# Replica Exchange log file") != 0) {
    mprinterr("Error: Expected '# Replica Exchange log file', got:\n%s\n", line.c_str());
    return 1;
  }

  // Read past metadata. Save expected number of exchanges.
  while (line[0] == '#') {
    line = buffer.GetLine();
    ++LineNum;
    if (line.empty()) {
      mprinterr("Error: No exchanges in rem log.\n");
      return 1;
    }
    ArgList columns( line );
    if (columns.hasKey("exchange")) break;
    mprintf("\t%s", line.c_str());
    if (columns.hasKey("numexchg")) {
      numexchg = columns.getNextInteger(-1);
    }
    if (columns.hasKey("Rep#,")) {
      if (columns[2] == "Neibr#,") type = HREMD;
      else if (columns[2] == "Velocity") type = TREMD;
    }
  }
  mprintf("\tRem log %s should contain %i exchanges\n", fname.c_str(), numexchg);
  if (numexchg < 1) {
    mprinterr("Error: Invalid number of exchanges (%i) in rem log.\n");
    return 1;
  }
  mprintf("\tExchange type is %s\n", ExchgDescription[type]);
  if (type == UNKNOWN) {
    mprinterr("Error: Could not identify exchange type.\n");
    return 1;
  }
  mprintf("\tFirst exchange is on line %i\n", LineNum);

  // Should currently be positioned at the first exchange. Need to read this
  // to determine how many replicas there are (and temperature for T-REMD).
  DataSet_RemLog::TmapType TemperatureMap;
  double t0;
  int n_replicas = 0;
  const char* ptr = buffer.Line();
  while (ptr != 0 && ptr[0] != '#') {
    if (type == TREMD) {
      // For temperature remlog create temperature map. Map temperatures to 
      // index + 1 since indices in the remlog start from 1.
      mprintf("DEBUG: Temp0= %s", ptr+32);
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
  mprintf("\t%i replicas.\n", n_replicas);
  if (n_replicas < 1) {
    mprinterr("Error: Detected less than 1 replica in remlog.\n");
    return 1;
  }
  DataSet* ds = datasetlist.AddSet( DataSet::REMLOG, fname, "remlog" );
  if (ds == 0) return 1;
  DataSet_RemLog& ensemble = static_cast<DataSet_RemLog&>( *ds );
  ensemble.AllocateReplicas(n_replicas);
  std::vector<DataSet_RemLog::ReplicaFrame> replicaFrames;
  std::vector<int> coordinateIndices;
  if (type == HREMD) {
    replicaFrames.resize( n_replicas );
    // All coord indices start equal to replica indices.
    // Indices start from 1 in remlogs (H-REMD only).
    coordinateIndices.resize( n_replicas );
    for (int replica = 0; replica < n_replicas; replica++)
      coordinateIndices[replica] = replica+1;
  } else if (type == TREMD)
    replicaFrames.resize(1);

  // Close and reopen the file, advance back to first exchange 
  buffer.CloseFile();
  buffer.OpenFileRead( fname );
  ptr = buffer.Line();
  while (ptr[0] == '#' && ptr[2] != 'e' && ptr[3] != 'x') ptr = buffer.Line();
  // Should now be positioned at 'exchange 1'.

  // Loop over all exchanges.
  ProgressBar progress( numexchg );
  for (int exchg = 0; exchg < numexchg; exchg++) {
    progress.Update( exchg );
    for (int replica = 0; replica < n_replicas; replica++) {
      // Read remlog line.
      ptr = buffer.Line();
      if (ptr == 0) {
        mprinterr("Error: reading remlog; unexpected end of file. Exchange=%i, replica=%i\n",
                  exchg+1, replica+1);
        return 1;
      }
      // ----- T-REMD ----------------------------
      if (type == TREMD) {
          if (replicaFrames[0].SetTremdFrame( ptr, TemperatureMap )) {
            mprinterr("Error reading TREMD line from rem log. Exchange=%i, replica=%i\n",
                      exchg+1, replica+1);
            return 1;
          }
          // Add replica frame to appropriate ensemble
          ensemble.Replica(replicaFrames[0].ReplicaIdx()-1).push_back( replicaFrames[0] );
      // ----- H-REMD ----------------------------
      } else if (type == HREMD) {
          if (replicaFrames[replica].SetHremdFrame( ptr, coordinateIndices[replica] )) {
            mprinterr("Error reading HREMD line from rem log. Exchange=%i, replica=%i\n",
                      exchg+1, replica+1);
            return 1;
          }
          // Add replica frame to appropriate ensemble
          ensemble.Replica(replica).push_back( replicaFrames[replica] );
      // -----------------------------------------
      } else {
        mprinterr("Error: remlog; unknown type.\n");
      }
    }
    if (type == HREMD) {
      // Determine whether exchanges occurred. Update coordinate indices accordingly.
      for (int replica = 0; replica < n_replicas; replica++) {
        if (replicaFrames[replica].Success()) {
          int partner = replicaFrames[replica].PartnerIdx() - 1;
          coordinateIndices[replica] = replicaFrames[partner].CoordsIdx();
        }
      }
    }
    // Read 'exchange N' line.
    ptr = buffer.Line();
  }

  buffer.CloseFile();
  // DEBUG - Print out replica 1 stats
  mprintf("Replica 1 Stats:\n"
          "%-10s %6s %6s %6s %12s %12s %12s S\n", "#Exchange", "RepIdx", "PrtIdx", "CrdIdx",
          "Temp0", "PE_X1", "PE_X2"); 
  for (DataSet_RemLog::replica_it it = ensemble.begin(0);
                                  it != ensemble.end(0); ++it)
    mprintf("%10u %6i %6i %6i %12.4f %12.4f %12.4f %1i\n", it - ensemble.begin(0) + 1,
            (*it).ReplicaIdx(), (*it).PartnerIdx(), (*it).CoordsIdx(), (*it).Temp0(), 
            (*it).PE_X1(), (*it).PE_X2(), (int)(*it).Success()); 
  return 0;
}
