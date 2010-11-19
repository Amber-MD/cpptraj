// TrajoutList
#ifndef INC_TRAJOUTLIST_H
#define INC_TRAJOUTLIST_H
//#include "TrajFile.h" // Frame.h
//#include "ArgList.h"
//#include "ParmFileList.h" // AmberParm.h
//#include "FrameList.h"
#include "CoordFileList.h"

class TrajoutList : public CoordFileList {

  public:
    
    TrajoutList();
    ~TrajoutList();
    // Add a traj file to the list with given access and associate with a parm
    // NOTE: worldsize is passed in as last arg to avoid include of PtrajMpi
    int Add(ArgList *A, ParmFileList *, int);
    int Write(int , Frame *, AmberParm*);
    void Close();
};
#endif

