#ifndef INC_COORDINATEINFO_H
#define INC_COORDINATEINFO_H
#include <string>
#include "ReplicaDimArray.h"
#include "Box.h"
#ifdef MPI
# include "Parallel.h"
#endif
/// All metadata associated with a Frame.
class CoordinateInfo {
  public:
    /// CONSTRUCTOR
    CoordinateInfo();
    /// CONSTRUCTOR - box, velocity, temperature, time
    CoordinateInfo(Box const&, bool, bool, bool);
    /// CONSTRUCTOR - box, coord, velocity, force, time TODO merge with above? 
    CoordinateInfo(Box const&, bool, bool, bool, bool);
    /// CONSTRUCTOR - ensemble size, remd dims, box, coords, velocity, force, temp., pH, redox, time, step, has repidx, has crdidx, use remd values
    CoordinateInfo(int, ReplicaDimArray const&, Box const&, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool);
    bool HasBox()              const { return box_.HasBox();            }
    const Box& TrajBox()       const { return box_;                     }
    int EnsembleSize()         const { return ensembleSize_;            }
    bool HasCrd()              const { return hasCrd_;                  }
    bool HasVel()              const { return hasVel_;                  }
    bool HasForce()            const { return hasFrc_;                  }
    bool HasTemp()             const { return hasTemp_;                 }
    bool Has_pH()              const { return has_pH_;                  }
    bool HasRedOx()            const { return hasRedox_;                }
    bool HasTime()             const { return hasTime_;                 }
    bool HasStep()             const { return hasStep_;                 }
    bool HasReplicaDims()      const { return (remdDim_.Ndims() != 0);  }
    bool HasRepIdx()           const { return hasrepidx_;               }
    bool HasCrdIdx()           const { return hascrdidx_;               }
    ReplicaDimArray const& ReplicaDimensions() const { return remdDim_; }
    bool UseRemdValues()       const { return useRemdValues_;           }
    void SetTime(bool m)        { hasTime_ = m; }
    void SetStep(bool s)        { hasStep_ = s; }
    void SetTemperature(bool t) { hasTemp_ = t; }
    void SetCrd(bool c)         { hasCrd_ = c;  }
    void SetVelocity(bool v)    { hasVel_ = v;  }
    void SetForce(bool f)       { hasFrc_ = f;  }
    void SetEnsembleSize(int s) { ensembleSize_ = s; }
    void SetBox(Box const& b)   { box_ = b;     }
    void SetReplicaDims(ReplicaDimArray const& r) { remdDim_ = r; }
    /// Print coordinate info to STDOUT
    void PrintCoordInfo(const char*, const char*) const;
    /// \return string containing info on present metadata
    std::string InfoString() const;
#   ifdef MPI
    int SyncCoordInfo(Parallel::Comm const&);
#   endif
    /// \return True if Frame would need to be re-setup based on CoordinateInfo
    bool operator !=(CoordinateInfo const& rhs) const {
      return (hasVel_ != rhs.hasVel_ ||
              hasFrc_ != rhs.hasFrc_ ||
              remdDim_.Ndims() != rhs.remdDim_.Ndims());
    }
  private:
    ReplicaDimArray remdDim_; ///< Hold info on any replica dimensions.
    Box box_;                 ///< Hold box information.
    int ensembleSize_;        ///< If coordinate ensemble, total # of replicas.
    bool hasCrd_;             ///< True if coords present. Now only relevant for NetCDF traj.
    bool hasVel_;             ///< True if coords have associated velocities.
    bool hasFrc_;             ///< True if coords have associated forces.
    bool hasTemp_;            ///< True if coords include temp info.
    bool has_pH_;             ///< True if coords include pH info.
    bool hasRedox_;           ///< True if coords include RedOx potential info.
    bool hasTime_;            ///< True if coords include time info.
    bool hasStep_;            ///< True if coords include step info.
    bool hasrepidx_;          ///< True if coords have replica indices
    bool hascrdidx_;          ///< True if coords have coordinate indices
    bool useRemdValues_;      ///< True if using remd_values in netcdf file
};
#endif
