#ifndef INC_TRAJIN_MULTI_H
#define INC_TRAJIN_MULTI_H
#include "Trajin.h"
#include "TrajIOarray.h"
#include "ReplicaInfo.h"
#include "Frame.h" // For RemdIdxType
// Forward declares
class ArgList;
class Topology;
class FileName;
/// Read in 1 frame at a time from multiple files.
class Trajin_Multi : public Trajin {
  public:
    Trajin_Multi();
    ~Trajin_Multi() { EndTraj(); }
    // ----- Inherited functions -------------------
    /// Set up trajectory for reading.
    int SetupTrajRead(FileName const&, ArgList&, Topology*);
    /// Read specified frame #.
    int ReadTrajFrame(int, Frame&);
    /// Prepare trajectory for reading.
    int BeginTraj();
    /// Close trajectory.
    void EndTraj();
    /// Print trajectory information.
    void PrintInfo(int) const;
    /// \return trajectory metadata.
    CoordinateInfo const& TrajCoordInfo() const { return cInfo_; }
    // ---------------------------------------------
  private:
    /// Type that will hold REMD indices
    typedef Frame::RemdIdxType RemdIdxType; // TODO put in ReplicaInfo

    TrajIOarray REMDtraj_;
    CoordinateInfo cInfo_; ///< Collective coord information for all replicas TODO Trajin?
    ReplicaInfo::TargetType targetType_; ///< Hold type of REMD frame being searched for.
    RemdIdxType remdtrajidx_; ///< Get frames with these indices on read
    double remdtrajtemp_;     ///< Get frames with this temperature on read
};
#endif 
