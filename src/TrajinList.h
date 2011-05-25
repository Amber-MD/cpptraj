// TrajinList
#ifndef INC_TRAJINLIST_H
#define INC_TRAJINLIST_H
#include "CoordFileList.h"

class TrajinList : public CoordFileList {

  public:
    
    TrajinList();
    ~TrajinList();
    // Add a traj file with given name and parm to the list
    int AddTrajin(char *, AmberParm *);
    // Add a traj file to the list based on input from arg list
    int Add(ArgList *, ParmFileList *);
    // TRAJIN: Set up frames to be processed 
    int SetupFrames();
};
#endif

