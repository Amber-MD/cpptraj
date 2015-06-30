#ifndef INC_TRAJ_GRO_H
#define INC_TRAJ_GRO_H
#include "TrajectoryIO.h"
#include "BufferedLine.h"
class Traj_Gro : public TrajectoryIO {
  public:
    Traj_Gro() : debug_(0), natom_(0), nbox_(0) {}
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_Gro(); }
    static void WriteHelp();
  private:
    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(std::string const&, Topology*);
    int setupTrajout(std::string const&, Topology*, CoordinateInfo const&,int, bool) { return 1;}
    int openTrajin() { return 1; }
    void closeTraj() {}
    int readFrame(int,Frame&) { return 1; }
    int writeFrame(int,Frame const&) { return 1; }
    void Info();
    int processWriteArgs(ArgList&) { return 0; }
    int readVelocity(int, Frame&) { return 1; }
    int processReadArgs(ArgList&)  { return 0; }

    double GetTimeValue(const char*) const;

    BufferedLine file_;
    int debug_;
    int natom_;
    int nbox_;
};
#endif
