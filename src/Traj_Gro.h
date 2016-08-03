#ifndef INC_TRAJ_GRO_H
#define INC_TRAJ_GRO_H
#include "TrajectoryIO.h"
#include "BufferedLine.h"
class Traj_Gro : public TrajectoryIO {
  public:
    Traj_Gro() : debug_(0), natom_(0), currentSet_(0), linesToRead_(0) {}
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_Gro(); }
    static void WriteHelp();
  private:
    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(FileName const&, Topology*);
    int setupTrajout(FileName const&, Topology*, CoordinateInfo const&,int, bool) { return 1;}
    int openTrajin();
    void closeTraj() { file_.CloseFile(); }
    int readFrame(int,Frame&);
    int writeFrame(int,Frame const&) { return 1; }
    void Info();
    int processWriteArgs(ArgList&) { return 0; }
    int readVelocity(int, Frame&)  { return 1; }
    int readForce(int, Frame&)     { return 1; }
    int processReadArgs(ArgList&)  { return 0; }

    double GetTimeValue(const char*) const;
    Box GetBox(const char*) const;

    BufferedLine file_;
    int debug_;
    int natom_;
    int currentSet_;
    int linesToRead_; ///< For blank reads
    FileName fname_; ///< File name TODO file_ should save this
    static const double GMX_VEL_TO_AMBER;
};
#endif
