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

    int AddRefFrame(Frame *, char *, AmberParm *,int);
    int AddFrame(Frame *, AmberParm *);
    AmberParm *GetFrameParm(int);
    int GetFrameIndex(char *);
    Frame *GetFrame(int idx);
    int ReplaceFrame(int, Frame *, AmberParm *);
    void Info();
    const char *FrameName(int);
    int NumFrames() { return Nframe; }
};
#endif
