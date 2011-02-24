#ifndef INC_FRAMELIST_H
#define INC_FRAMELIST_H
#include "Frame.h"
#include "ParmFileList.h"
#include <vector>
#include <string>
// FrameList
class FrameList {
    std::vector<Frame*> frameList;
    ParmFileList FrameParm;
    std::vector<std::string> frameNames;
    std::vector<int> frameNums;
    int Nframe;
  
  public:

    FrameList();
    ~FrameList();

    int Add(Frame *, char *, AmberParm *,int);
    AmberParm *GetFrameParm(int);
    int GetFrameIndex(char *);
    Frame *GetFrame(int idx);
    void Info();
};
#endif
