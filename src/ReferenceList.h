// ReferenceList
#ifndef INC_REFERENCELIST_H
#define INC_REFERENCELIST_H
#include "CoordFileList.h"
#include "FrameList.h"
class ReferenceList : public CoordFileList {
    std::vector<bool> Average;
  public:
    
    ReferenceList();
    ~ReferenceList();
    // Add a traj file to the list with given access and associate with a parm
    int Add(char*,ArgList *A, AmberParm *);
    // REFERENCE: Set up frames to be processed
    int SetupRefFrames(FrameList *);
};
#endif

