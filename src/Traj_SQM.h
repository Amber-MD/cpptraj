#ifndef INC_TRAJ_SQM_H
#define INC_TRAJ_SQM_H
#include "TrajectoryIO.h"
/// Write out sqm input file.
class Traj_SQM : public TrajectoryIO {
  public:
    Traj_SQM() : singleWrite_(false), chargeIsSet_(false), charge_(0), sqmParm_(0) {}
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_SQM(); }
    static void WriteHelp();
  private:
    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&)               { return false; }
    int processReadArgs(ArgList&)                  { return 0;     }
    int setupTrajin(FileName const&, Topology*)    { return 1;     }
    int openTrajin()                               { return 1;     }
    int readFrame(int,Frame&)                      { return 1;     }
    int readVelocity(int, Frame&)                  { return 1;     }
    int readForce(int, Frame&)                     { return 1;     }
    void closeTraj()                               { return;       }
    int processWriteArgs(ArgList&, DataSetList const&);
    int setupTrajout(FileName const&, Topology*, CoordinateInfo const&,int, bool);
    int writeFrame(int,Frame const&);
    void Info();

    bool singleWrite_;
    bool chargeIsSet_;
    int charge_;
    Topology* sqmParm_;
    CpptrajFile outfile_;
    std::string header_;
};
#endif
