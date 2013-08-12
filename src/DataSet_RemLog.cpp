#include <cstdio> // sscanf
#include "DataSet_RemLog.h"
#include "CpptrajStdio.h"

DataSet_RemLog::DataSet_RemLog() :
  DataSet(REMLOG, 10, 4, 0) // 0 dim indicates DataSet-specific write 
{}

void DataSet_RemLog::AllocateReplicas(int n_replicas) {
  ensemble_.clear();
  ensemble_.resize( n_replicas );
}

int DataSet_RemLog::NumExchange() const {
  if (ensemble_.empty())
    return 0;
  else // Each member of the ensemble should have same # exchanges.
    return (int)ensemble_[0].size();
}

bool DataSet_RemLog::ValidEnsemble() const {
  ReplicaEnsemble::const_iterator member = ensemble_.begin();
  size_t first_size = (*member).size();
  for (; member != ensemble_.end(); ++member) {
    if ((*member).size() != first_size) {
      mprinterr("Error: In remlog data set %s size of ensemble member %zu (%zu) !="
                " size of first member (%zu)\n", Name().c_str(), // TODO: Change to legend
                member - ensemble_.begin() + 1, (*member).size(), first_size);
      return false;
    }
  }
  return true;
}

void DataSet_RemLog::TrimLastExchange() {
  if (ensemble_.empty()) return;
  ReplicaEnsemble::iterator member = ensemble_.begin();
  size_t min_size = (*member).size();
  ++member;
  for (; member != ensemble_.end(); ++member) {
    if ((*member).size() < min_size) min_size = (*member).size();
  }
  // Resize all member arrays to minimum
  for (member = ensemble_.begin(); member != ensemble_.end(); ++member)
    (*member).resize( min_size );
}

// -----------------------------------------------------------------------------
/* Format:
 * '(i2,6f10.2,i8)'
# Rep#, Velocity Scaling, T, Eptot, Temp0, NewTemp0, Success rate (i,i+1), ResStruct#
 1     -1.00      0.00   -433.24    300.00    300.00      0.00      -1
 * Order during REMD is exchange -> MD, so NewTemp0 is the temp. that gets
 * simulated.
 */
int DataSet_RemLog::ReplicaFrame::SetTremdFrame(const char* ptr, TmapType const& TemperatureMap)
{
  double scaling;
  if ( sscanf(ptr, "%2i%10lf%*10f%10lf%*10f%10lf", &coordsIdx_, &scaling, &PE_x1_,
              &temp0_) != 4 )
    return 1;
  // Figure out replica index from temperature map.
  TmapType::const_iterator tmap = TemperatureMap.find( temp0_ );
  if (tmap == TemperatureMap.end()) {
    mprinterr("Error: replica temperature %.2f not found in temperature map.\n", temp0_);
    return 1;
  }
  replicaIdx_ = (*tmap).second;
  // Partner index is not used in T-REMD.
  partnerIdx_ = 0; 
  // Exchange occurred if velocity scaling is > 0
  success_ = (scaling > 0.0);

  return 0;
}

/* Format:
 * '(2i6,5f10.2,4x,a,2x,f10.2)'
     1     8    300.00 -25011.03 -24959.58    -27.48      0.00    F        0.00
 */
int DataSet_RemLog::ReplicaFrame::SetHremdFrame( const char* ptr, int currentCoordIdx ) {
  if ( sscanf(ptr, "%6i%6i%10lf%10lf%10lf", &replicaIdx_, &partnerIdx_,
              &temp0_, &PE_x1_, &PE_x2_) != 5 )
    return 1;
  coordsIdx_ = currentCoordIdx;
  // Check if an exchange occured this frame.
  switch ( ptr[66] ) {
    case 'T': success_ = true;  break;
    case 'F': success_ = false; break;
    default: // Should only get here with malformed HREMD log file.
      mprinterr("Error: expected only 'T' or 'F' at character 67, got %c\n", ptr[66]);
      return 1;
  }
  return 0;
}
