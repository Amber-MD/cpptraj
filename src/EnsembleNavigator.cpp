#include "EnsembleNavigator.h"
#include "CpptrajStdio.h"

int EnsembleNavigator::AddEnsembles(TrajinList::ensemble_it const& begin,
                                    TrajinList::ensemble_it const& end)
{
  FirstParm_ = 0;
  for (TrajinList::ensemble_it ens = begin; ens != end; ++ens)
  {
    if (FirstParm_ == 0) {
      FirstParm_ = (*ens)->Traj().Parm();
      ensCoordInfo_ = (*ens)->EnsembleCoordInfo();
    } else {
      if (FirstParm_ != (*ens)->Traj().Parm()) {
        mprinterr("Error: Trajectory parallelized 'ensemble' currently requires all\n"
                  "Error:   ensembles use the same topology file.\n");
        return 1;
      }
      // Want to ensure frame is properly allocated for all ensembles, so
      // e.g. if one ensemble has velocity info allocate for all.
      if (ensCoordInfo_.HasVel() != (*ens)->EnsembleCoordInfo().HasVel())
        ensCoordInfo_.SetVelocity( true );
      if (ensCoordInfo_.HasForce() != (*ens)->EnsembleCoordInfo().HasForce())
        ensCoordInfo_.SetForce( true );
      // Sanity check
      if (ensCoordInfo_.ReplicaDimensions().Ndims() !=
          (*ens)->EnsembleCoordInfo().ReplicaDimensions().Ndims())
      {
        mprinterr("Internal Error: Replica dimensions changed.\n");
        return 1;
      }
    }
    IDX_.AddTraj( (*ens)->Traj().Counter().TotalReadFrames(),
                  (*ens)->Traj().Counter().Start(),
                  (*ens)->Traj().Counter().Offset() );
    Ensembles_.push_back( *ens );
  }
  return 0;
}
