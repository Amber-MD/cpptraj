#ifndef INC_TRAJ_CHARMMCOR_H
#define INC_TRAJ_CHARMMCOR_H
#include "TrajectoryIO.h"
#include "CpptrajFile.h"
/// Read CHARMM Cor file
class Traj_CharmmCor : public TrajectoryIO {
  public:
    Traj_CharmmCor() : corAtom_(0), extendedFmt_(false) {}
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_CharmmCor(); }
    //static void ReadHelp();
  private:
    // TrajectoryIO functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(FileName const&, Topology*);
    int openTrajin();
    void closeTraj();
    int readFrame(int,Frame&);
    int setupTrajout(FileName const&, Topology*, CoordinateInfo const&,int, bool);
    int writeFrame(int,Frame const&);
    void Info();
    int processWriteArgs(ArgList&, DataSetList const&) { return 0; }
    int readVelocity(int, Frame&)  { return 1; }
    int readForce(int, Frame&)     { return 1; }
    int processReadArgs(ArgList&)  { return 0; }

    CpptrajFile file_;
    int corAtom_;      ///< # of atoms in Cor file.
    bool extendedFmt_; ///< True for wide columns
};
#endif
