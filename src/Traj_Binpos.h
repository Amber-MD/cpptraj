#ifndef INC_TRAJ_BINPOS
#define INC_TRAJ_BINPOS
#include "TrajectoryIO.h"
// Class: Traj_Binpos
/// Read and write SCRIPPS BINPOS format
class Traj_Binpos : public TrajectoryIO {
  public:
    Traj_Binpos();
    ~Traj_Binpos();
  private:
    // Inherited functions
    bool ID_TrajFormat();
    int setupTrajin(Topology*);
    int setupTrajout(Topology*,int);
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    void info();
    //int processWriteArgs(ArgList*);

    int bpnatom_;
    int bpnatom3_;
    size_t frameSize_;
    float* bpbuffer_;
};
#endif
