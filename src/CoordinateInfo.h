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
    CoordinateInfo() : ensembleSize_(0), hasCrd_(true), hasVel_(false),
                       hasTemp_(false), hasTime_(false), hasFrc_(false) {}
    /// CONSTRUCTOR - box, velocity, temperature, time
    CoordinateInfo(Box const& b, bool v, bool t, bool m) :
      box_(b), ensembleSize_(0), hasCrd_(true), hasVel_(v),
      hasTemp_(t), hasTime_(m), hasFrc_(false) {}
    /// CONSTRUCTOR - all except ensemble size
    CoordinateInfo(ReplicaDimArray const& r, Box const& b, bool v, bool t, bool m, bool f) :
      remdDim_(r), box_(b), ensembleSize_(0), hasCrd_(true), hasVel_(v),
      hasTemp_(t), hasTime_(m), hasFrc_(f) {}
    /// CONSTRUCTOR - All
    CoordinateInfo(int e, ReplicaDimArray const& r, Box const& b, bool v, bool t, bool m, bool f) :
      remdDim_(r), box_(b), ensembleSize_(e), hasCrd_(true), hasVel_(v),
      hasTemp_(t), hasTime_(m), hasFrc_(f) {}
    bool HasBox()              const { return box_.HasBox();            }
    const Box& TrajBox()       const { return box_;                     }
    int EnsembleSize()         const { return ensembleSize_;            }
    bool HasCrd()              const { return hasCrd_;                  }
    bool HasVel()              const { return hasVel_;                  }
    bool HasTemp()             const { return hasTemp_;                 }
    bool HasTime()             const { return hasTime_;                 }
    bool HasForce()            const { return hasFrc_;                  }
    bool HasReplicaDims()      const { return (remdDim_.Ndims() != 0);  }
    ReplicaDimArray const& ReplicaDimensions() const { return remdDim_; }
    void SetTime(bool m)        { hasTime_ = m; }
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
    bool hasTemp_;            ///< True if coords include temp info.
    bool hasTime_;            ///< True if coords include time info.
    bool hasFrc_;             ///< True if coords have associated forces.
};
#endif
