#include "EnsembleIn.h"
#include "CpptrajStdio.h"
#include "Constants.h" // FNE

EnsembleIn::EnsembleIn() :
  targetType_(ReplicaInfo::NONE), badEnsemble_(0), debug_(0)
# ifdef MPI
  , member_(-1)
# endif
{}

#ifdef MPI
int EnsembleIn::GatherTemperatures(double* tAddress, std::vector<double>& allTemps,
                                   Parallel::Comm const& commIn)
{
  if (commIn.AllGather(tAddress, 1, MPI_DOUBLE, &allTemps[0])) {
    rprinterr("Error: Gathering temperatures.\n");
    return 1; // TODO: Better parallel error check
  }
  return 0;
}

int EnsembleIn::GatherIndices(int* iAddress, std::vector<RemdIdxType>& allIndices,
                              int Ndims, Parallel::Comm const& commIn)
{
  std::vector<int> all_indices(allIndices.size() * Ndims, 0);
  if (commIn.AllGather(iAddress, Ndims, MPI_INT, &all_indices[0])) {
    rprinterr("Error: Gathering replica indices\n");
    return 1; // TODO: Better parallel error check
  }
  std::vector<int>::const_iterator idx_it = all_indices.begin();
  for (std::vector<RemdIdxType>::iterator it = allIndices.begin();
                                          it != allIndices.end();
                                        ++it, idx_it += Ndims)
    it->assign(idx_it, idx_it + Ndims);
  return 0;
}

#ifdef TIMER
double EnsembleIn::total_mpi_allgather_ = 0.0;
double EnsembleIn::total_mpi_sendrecv_ = 0.0;

void EnsembleIn::TimingData(double trajin_time) {
  if (total_mpi_allgather_ > 0.0 || total_mpi_sendrecv_ > 0.0) {
    double other_time = trajin_time - total_mpi_allgather_ - total_mpi_sendrecv_;
    rprintf("MPI_TIME:\tallgather: %.4f s (%.2f%%), sendrecv: %.4f s (%.2f%%), Other:  %.4f s (%.2f%%)\n",
            total_mpi_allgather_, (total_mpi_allgather_ / trajin_time)*100.0,
            total_mpi_sendrecv_,  (total_mpi_sendrecv_  / trajin_time)*100.0,
            other_time, (other_time / trajin_time)*100.0 );
  }
}
#endif /* TIMER */
#endif /* MPI */

void EnsembleIn::PrintReplicaInfo() const {
  if (targetType_ == ReplicaInfo::TEMP) {
    mprintf("  Ensemble Temperature Map:\n");
    for (ReplicaInfo::Map<double>::const_iterator tmap = TemperatureMap_.begin();
                                            tmap != TemperatureMap_.end(); ++tmap)
      mprintf("\t%10.2f -> %i\n", tmap->first, tmap->second);
  } else if (targetType_ == ReplicaInfo::INDICES) {
    mprintf("  Ensemble Indices Map:\n");
    for (ReplicaInfo::Map<Frame::RemdIdxType>::const_iterator imap = IndicesMap_.begin();
                                                        imap != IndicesMap_.end(); ++imap)
    {
      mprintf("\t{");
      for (Frame::RemdIdxType::const_iterator idx = imap->first.begin();
                                              idx != imap->first.end(); ++idx)
        mprintf(" %i", *idx);
      mprintf(" } -> %i\n", imap->second);
    }
  }
}

int EnsembleIn::SetTemperatureMap(std::vector<double> const& allTemps) {
  if (TemperatureMap_.CreateMap( allTemps )) {
    rprinterr("Error: Ensemble: Duplicate temperature detected (%.2f) in ensemble %s\n"
              "Error:   If this is an H-REMD ensemble try the 'nosort' keyword.\n",
              TemperatureMap_.Duplicate(), traj_.Filename().full());
    return 1;
  }
  return 0;
}

int EnsembleIn::SetIndicesMap(std::vector<RemdIdxType> const& allIndices) {
  if (IndicesMap_.CreateMap( allIndices )) {
    rprinterr("Error: Ensemble: Duplicate indices detected in ensemble %s:",
              traj_.Filename().full());
    for (RemdIdxType::const_iterator idx = IndicesMap_.Duplicate().begin();
                                     idx != IndicesMap_.Duplicate().end(); ++idx)
      rprinterr(" %i", *idx);
    rprinterr("\n");
    return 1;
  }
  return 0;
}

int EnsembleIn::SetIdxValMap(std::vector<RemdIdxType> const& allIndices,
                             std::vector<Darray> const& allValues)
{
  if (allIndices.size() != allValues.size()) {
    rprinterr("Internal Error: EnsembleIn::SetIdxValMap(): allIndices size %zu != allValues size %zu\n",
              allIndices.size(), allValues.size());
    return 1;
  }

  if (debug_ > 0) {
    mprintf("DEBUG: index/value map:\n");
    for (unsigned int member = 0; member != allIndices.size(); member++) {
      mprintf("\t");
      for (Frame::RemdIdxType::const_iterator idx = allIndices[member].begin();
                                              idx != allIndices[member].end(); ++idx)
        mprintf(" %6i", *idx);
      mprintf(" : ");
      for (Darray::const_iterator dval = allValues[member].begin();
                                  dval != allValues[member].end(); ++dval)
        mprintf(" %8.3f", *dval);
      mprintf("\n");
    }
  }

  typedef std::pair<RemdIdxType, Darray> ivPair;
  for (unsigned int idx = 0; idx != allIndices.size(); idx++)
    RemdIdxValMap_.insert( ivPair(allIndices[idx], allValues[idx]) );
  return 0;
}

static inline void indexprint(Frame::RemdIdxType const& IDXin) {
  for (Frame::RemdIdxType::const_iterator idx = IDXin.begin();
                                          idx != IDXin.end(); ++idx)
    mprintf(" %6i", *idx);
}

static inline void valprint(std::vector<double> const& VALin) {
  for (std::vector<double>::const_iterator dval = VALin.begin();
                                           dval != VALin.end(); ++dval)
    mprintf(" %8.3f", *dval);
}

/// Floating point not equals.
static inline bool FNE(double v1, double v2) {
  double delta = v1 - v2;
  if (delta < 0.0) delta = -delta;
  return (delta > Constants::SMALL);
}

// Compare Darrays
static inline bool DNE(std::vector<double> const& lhs,
                       std::vector<double> const& rhs)
{
  if (lhs.size() != rhs.size()) return true;
  for (unsigned int idx = 0; idx < lhs.size(); idx++)
    if (FNE(lhs[idx], rhs[idx])) return true;
  return false;
}

int EnsembleIn::CheckIdxValMap(RemdIdxValMapType const& mapIn) const {
  if (mapIn.size() != RemdIdxValMap_.size()) {
    mprintf("Warning: REMD index/value map size %zu does not match previous map size %zu\n",
            RemdIdxValMap_.size(), mapIn.size());
    return 1;
  }
  int match = 0;
  for (RemdIdxValMapType::const_iterator it = RemdIdxValMap_.begin();
                                         it != RemdIdxValMap_.end(); ++it)
  {
    RemdIdxValMapType::const_iterator jt = mapIn.find( it->first );
    if (jt == mapIn.end()) {
      mprintf("Warning: Index");
      indexprint( it->first );
      mprintf(" for this ensemble not found in previous ensemble.\n");
      match = 1;
    } else if (DNE(jt->second, it->second)) {
      mprintf("Warning: Index");
      indexprint( it->first );
      mprintf(" values ");
      valprint( it->second );
      mprintf(" do not match previous map values ");
      valprint( jt->second );
      mprintf("\n");
      match = 1;
    }
  }
  return match;
}
