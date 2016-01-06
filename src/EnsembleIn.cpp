#include "EnsembleIn.h"
#include "CpptrajStdio.h"

#ifdef MPI
int EnsembleIn::SetupComms(int ngroups) {
  // Make sure that ngroups is a multiple of total # threads.
  if ( (Parallel::World().Size() % ngroups) != 0 ) {
    mprinterr("Error: # of threads (%i) must be a multiple of # replicas (%i)\n",
              Parallel::World().Size(), ngroups);
    return 1;
  }
  // Split into comms for parallel across trajectory
  int ID = Parallel::World().Rank() / (Parallel::World().Size() / ngroups);
  rprintf("DEBUG: My trajComm ID is %i\n", ID);
  trajComm_ = Parallel::World().Split( ID );
  // Ensemble comm is orthogonal to trajectory comm
  ID = Parallel::World().Rank() % trajComm_.Size();
  rprintf("DEBUG: My ensembleComm ID is %i\n", ID);
  ensembleComm_ = Parallel::World().Split( ID );
  // Create comm from masters of trajComm_ for setting up ensemble.
  if (trajComm_.Master())
    ID = 0;
  else
    ID = MPI_UNDEFINED;
  if (ID == 0)
    rprintf("DEBUG: I am in the ensemble setup comm\n");
  else
    rprintf("DEBUG: Not part of ensemble setup comm\n");
  setupComm_ = Parallel::World().Split( ID );
  return 0;
}

int EnsembleIn::GatherTemperatures(double* tAddress, std::vector<double>& allTemps) const {
  if (setupComm_.AllGather(tAddress, 1, MPI_DOUBLE, &allTemps[0])) {
    rprinterr("Error: Gathering temperatures.\n");
    return 1; // TODO: Better parallel error check
  }
  return 0;
}

int EnsembleIn::GatherIndices(int* iAddress, std::vector<RemdIdxType>& allIndices, int Ndims) const
{
  std::vector<int> all_indices(allIndices.size() * Ndims, 0);
  if (setupComm_.AllGather(iAddress, Ndims, MPI_INT, &all_indices[0])) {
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
#endif
#endif /* MPI */

void EnsembleIn::PrintReplicaInfo() const {
  if (targetType_ == ReplicaInfo::TEMP) {
    mprintf("  Ensemble Temperature Map:\n");
    for (ReplicaMap<double>::const_iterator tmap = TemperatureMap_.begin();
                                            tmap != TemperatureMap_.end(); ++tmap)
      mprintf("\t%10.2f -> %i\n", tmap->first, tmap->second);
  } else if (targetType_ == ReplicaInfo::INDICES) {
    mprintf("  Ensemble Indices Map:\n");
    for (ReplicaMap<Frame::RemdIdxType>::const_iterator imap = IndicesMap_.begin();
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
