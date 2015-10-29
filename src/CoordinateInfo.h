#ifndef INC_COORDINATEINFO_H
#define INC_COORDINATEINFO_H
#include "ReplicaDimArray.h"
#include "Box.h"
/// All metadata associated with a Frame.
class CoordinateInfo {
  public:
    /// CONSTRUCTOR
    CoordinateInfo() : ensembleSize_(0), hasVel_(false), hasTemp_(false), hasTime_(false), hasFrc_(false) {}
    /// CONSTRUCTOR - box, velocity, temperature, time
    CoordinateInfo(Box const& b, bool v, bool t, bool m) :
      box_(b), ensembleSize_(0), hasVel_(v), hasTemp_(t), hasTime_(m), hasFrc_(false) {}
    /// CONSTRUCTOR - all except ensemble size
    CoordinateInfo(ReplicaDimArray const& r, Box const& b, bool v, bool t, bool m, bool f) :
      remdDim_(r), box_(b), ensembleSize_(0), hasVel_(v), hasTemp_(t), hasTime_(m), hasFrc_(f) {}
    /// CONSTRUCTOR - All
    CoordinateInfo(int e, ReplicaDimArray const& r, Box const& b, bool v, bool t, bool m, bool f) :
      remdDim_(r), box_(b), ensembleSize_(e), hasVel_(v), hasTemp_(t), hasTime_(m), hasFrc_(f) {}
    bool HasBox()              const { return box_.HasBox();            }
    const Box& TrajBox()       const { return box_;                     }
    int EnsembleSize()         const { return ensembleSize_;            }
    bool HasVel()              const { return hasVel_;                  }
    bool HasTemp()             const { return hasTemp_;                 }
    bool HasTime()             const { return hasTime_;                 }
    bool HasForce()            const { return hasFrc_;                  }
    bool HasReplicaDims()      const { return (remdDim_.Ndims() != 0);  }
    ReplicaDimArray const& ReplicaDimensions() const { return remdDim_; }
    void SetTime(bool m)        { hasTime_ = m; }
    void SetTemperature(bool t) { hasTemp_ = t; }
    void SetVelocity(bool v)    { hasVel_ = v;  }
    void SetForce(bool f)       { hasFrc_ = f;  }
    void SetEnsembleSize(int s) { ensembleSize_ = s; }
    void SetBox(Box const& b)   { box_ = b;     }
    void SetReplicaDims(ReplicaDimArray const& r) { remdDim_ = r; }
  private:
    ReplicaDimArray remdDim_; ///< Hold info on any replica dimensions.
    Box box_;                 ///< Hold box information.
    int ensembleSize_;        ///< If coordinate ensemble, total # of replicas.
    bool hasVel_;             ///< True if coords have associated velocities.
    bool hasTemp_;            ///< True if coords include temp info.
    bool hasTime_;            ///< True if coords include time info.
    bool hasFrc_;             ///< True if coords have associated forces.
};
#endif
