#ifndef TRAJIN_SINGLE_H
#define TRAJIN_SINGLE_H
#include "Trajin.h"
/// Class for reading in single trajectories.
class Trajin_Single : public Trajin {
  public:
    Trajin_Single();
    ~Trajin_Single();

    int SetupTrajRead(std::string const&, ArgList *, Topology *, bool);

    int SetupTrajRead(std::string const&, ArgList *, Topology *);
    int BeginTraj(bool);
    void EndTraj();
    int GetNextFrame(Frame&);
    void PrintInfo(int);
    bool HasVelocity();
  private:
    TrajectoryIO* trajio_; ///< Hold class that will interface with traj format.
    TrajectoryIO* velio_;  ///< Hold class that will interface with opt. mdvel file.
    bool trajIsOpen_;      ///< True is trajectory is open. 
};
#endif
