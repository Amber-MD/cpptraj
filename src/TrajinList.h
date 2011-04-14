// TrajinList
#ifndef INC_TRAJINLIST_H
#define INC_TRAJINLIST_H
#include "CoordFileList.h"

class TrajinList : public CoordFileList {

  public:
    
    TrajinList();
    ~TrajinList();
    // Add a traj file to the list with given access and associate with a parm
    // NOTE: worldsize is passed in as last arg to avoid include of PtrajMpi
    int Add(ArgList *A, ParmFileList *, int);
    // TRAJIN: Set up frames to be processed 
    int SetupFrames(int,int);
};
#endif

