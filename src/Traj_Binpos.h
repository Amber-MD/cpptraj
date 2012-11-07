#ifndef INC_TRAJ_BINPOS
#define INC_TRAJ_BINPOS
#include "TrajectoryIO.h"
// Class: Traj_Binpos
/// Read and write SCRIPPS BINPOS format
class Traj_Binpos : public TrajectoryIO {
  public:
    Traj_Binpos();
    static TrajectoryIO* Alloc() { return (TrajectoryIO*)new Traj_Binpos(); }
    ~Traj_Binpos();
  private:
    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(std::string const&, Topology*);
    int setupTrajout(std::string const&, Topology*, int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    void Info();
    int processWriteArgs(ArgList&) { return 0; }

    int bpnatom_;
    int bpnatom3_;
    size_t frameSize_;
    float* bpbuffer_;
    CpptrajFile file_;

    int readVelocity(int, double*) { return 1; }
    int readIndices(int,int*) { return 1; }
    int processReadArgs(ArgList&) { return 0; }
};
#endif
