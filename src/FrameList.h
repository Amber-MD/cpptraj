#ifndef INC_FRAMELIST_H
#define INC_FRAMELIST_H
#include "Frame.h"
#include "ArgList.h"
#include "ParmFileList.h"

class FrameList {
    Frame **frameList;
    ParmFileList FrameParm;
    ArgList frameNames;
    int Nframe;
  
  public:

    FrameList();
    ~FrameList();

    int Add(Frame *, char *, AmberParm *);
    //Frame *GetFrame(char *);
    AmberParm *GetFrameParm(int);
    int GetFrameIndex(char *);
    Frame *GetFrame(int idx);
    void Info();
};
#endif
