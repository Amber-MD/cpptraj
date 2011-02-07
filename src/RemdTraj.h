#ifndef INC_REMDTRAJ_H
#define INC_REMDTRAJ_H
//#include "TrajFile.h"
//#include "TrajFileList.h"
#include "TrajinList.h"
//#include "ArgList.h"

class RemdTraj : public TrajFile {
    double remdtrajtemp;
    TrajinList REMDtraj;
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
