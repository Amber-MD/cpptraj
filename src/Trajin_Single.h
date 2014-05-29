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
    /// Set up trajectory for reading, check box info.
    int SetupTrajRead(std::string const&, ArgList&, Topology*);
    int BeginTraj(bool);
    void EndTraj();
    int ReadTrajFrame(int, Frame&);
    void PrintInfo(int) const;
    bool HasVelocity() const;
    /// \return Any replica dimension information present.
    ReplicaDimArray const& TrajReplicaDimInfo() const;
    int EnsembleSize() const { return 0; }
  private:
    TrajectoryIO* trajio_; ///< Hold class that will interface with traj format.
    TrajectoryIO* velio_;  ///< Hold class that will interface with opt. mdvel file.
    bool trajIsOpen_;      ///< True if trajectory is open.
    static const ReplicaDimArray emptyReplicaDimArray_; ///< When no replica dim info present 
};
#endif
