#ifndef INC_TRAJ_SQM_H
#define INC_TRAJ_SQM_H
#include "TrajectoryIO.h"
/// Write out sqm input file.
class Traj_SQM : public TrajectoryIO {
  public:
    Traj_SQM() : singleWrite_(false), sqmParm_(0) {}
    static TrajectoryIO* Alloc() { return (TrajectoryIO*)new Traj_SQM(); }
  private:
    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&)               { return false; }
    int processReadArgs(ArgList&)                  { return 0;     }
    int setupTrajin(std::string const&, Topology*) { return 1;     }
    int openTrajin()                               { return 1;     }
    int readFrame(int,Frame&)                      { return 1;     }
    int readVelocity(int, Frame&)                  { return 1;     }
    void closeTraj()                               { return;       }
    int processWriteArgs(ArgList&);
    int setupTrajout(std::string const&, Topology*, int, bool);
    int writeFrame(int,Frame const&);
    void Info();

    bool singleWrite_;
    Topology* sqmParm_;
    CpptrajFile outfile_;
};
#endif
