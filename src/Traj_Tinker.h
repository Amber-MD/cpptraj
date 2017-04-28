#ifndef INC_TRAJ_TINKER_H
#define INC_TRAJ_TINKER_H
#include "TrajectoryIO.h"
#include "TinkerFile.h"
/// TrajectoryIO class for reading coordinates from Tinker XYZ/ARC files.
class Traj_Tinker : public TrajectoryIO {
  public:
    Traj_Tinker();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_Tinker(); }
  private:
    int currentSet_;
    TinkerFile file_; 

    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(FileName const&, Topology*);
    int setupTrajout(FileName const&, Topology*, CoordinateInfo const&,int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,Frame&);
    int writeFrame(int,Frame const&);
    void Info();
    int processWriteArgs(ArgList&) { return 0; }
    int readVelocity(int, Frame&)  { return 1; }
    int readForce(int, Frame&)     { return 1; }
    int processReadArgs(ArgList&)  { return 0; }
};
#endif
