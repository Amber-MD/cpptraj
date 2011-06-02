#ifndef INC_REMDTRAJ_H
#define INC_REMDTRAJ_H
#include "TrajinList.h"
#include "TrajoutList.h"

class RemdTraj : public TrajFile {
    double remdtrajtemp;
    TrajinList REMDtraj;
    ArgList *RemdOutArgs;
    TrajoutList REMDtrajout;
    double *TemperatureList;
    Frame *remdframe;
    char *replicaName;
    int numReplicas;
  public:
    RemdTraj();
    ~RemdTraj();

    void SetReplicaName(char*, ArgList*);
    int NoTempInfo(TrajFile *);

    int SetupRead();
    int open();
    void close();
    int getFrame(int);
    void Info();
};
#endif
