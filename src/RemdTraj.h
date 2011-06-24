#ifndef INC_REMDTRAJ_H
#define INC_REMDTRAJ_H
/// Class: RemdTraj
/// Used to read in replica trajectories. Hold each replica trajectory
/// in a CoordFileList. During reads, only the frame with a temperature
/// matching remdtrajtemp will be used.
#include "TrajinList.h"
#include "TrajoutList.h"

class RemdTraj : public TrajectoryFileBase {
    double remdtrajtemp;
    TrajinList REMDtraj;
    ArgList *RemdOutArgs;
    TrajoutList REMDtrajout;
    double *TemperatureList;
    Frame *remdframe;
    int numReplicas;

    bool NoTempInfo(TrajectoryFile *);
  public:
    RemdTraj();
    ~RemdTraj();
    // Inherited functions
    int SetupRead(char*,ArgList*,AmberParm*);
    int SetupFrameInfo();
    int BeginTraj(bool);
    int EndTraj();
    int GetNextFrame(double*, double*,double*);
    void PrintInfo(int);
};
#endif
