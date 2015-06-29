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
    int setupTrajout(std::string const&, Topology*, CoordinateInfo const&,int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,Frame&);
    int writeFrame(int,Frame const&);
    void Info();
    int processWriteArgs(ArgList&);
    int readVelocity(int, Frame&);
    int processReadArgs(ArgList&)  { return 0; }

    double GetTimeValue(const char*) const;

    BufferedLine file_;
    int debug_;
    int natom_;
    int nbox_;
};
#endif
