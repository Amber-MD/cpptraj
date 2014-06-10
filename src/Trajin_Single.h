#ifndef TRAJIN_SINGLE_H
#define TRAJIN_SINGLE_H
#include "Trajin.h"
/// Class for reading in single trajectories.
class Trajin_Single : public Trajin {
  public:
    Trajin_Single();
    ~Trajin_Single();
    /// Set up trajectory for reading, optionally checking box info.
    int SetupTrajRead(std::string const&, ArgList&, Topology*, bool);
    // ----- Inherited functions -----------------
    /// Set up trajectory for reading, check box info.
    int SetupTrajRead(std::string const&, ArgList&, Topology*);
    int ReadTrajFrame(int, Frame&);
    int BeginTraj(bool);
    void EndTraj();
    void PrintInfo(int) const;
    bool HasVelocity() const;
    /// \return Any replica dimension information present.
    ReplicaDimArray const& TrajReplicaDimInfo() const {return trajRepDimInfo_;}
    int EnsembleSize() const { return 0; }
    // NOTE: The following are currently for testing Trajin_Ensemble
    void EnsembleInfo() const {} 
    int EnsembleSetup(FrameArray&,FramePtrArray&) {return 1;}
    int ReadEnsemble(int,FrameArray&,FramePtrArray&) {return 1;}
    bool  BadEnsemble() const { return true; }
  private:
    TrajectoryIO* trajio_; ///< Hold class that will interface with traj format.
    TrajectoryIO* velio_;  ///< Hold class that will interface with opt. mdvel file.
    bool trajIsOpen_;      ///< True if trajectory is open.
    ReplicaDimArray trajRepDimInfo_;
};
#endif
