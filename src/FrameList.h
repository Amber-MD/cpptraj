#ifndef INC_FRAMELIST_H
#define INC_FRAMELIST_H
#include <vector>
#include <string>
#include "Frame.h"
#include "ParmFileList.h"
// Class: FrameList
/// Hold a series of Frame classes. 
/** Optionally hold a corresponding name and number (in the case of 
  * reference frames, this is the name of the reference coordinate 
  * file and the frame number). 
  */
// NOTE: Eventually store a vector of Frames, not Frame*s
class FrameList {
    std::vector<Frame*> frameList;
    ParmFileList FrameParm;
    std::vector<std::string> frameNames;
    std::vector<std::string> frameTags;
    std::vector<int> frameNums;
    int Nframe;
    int referenceFrameNum;
  
  public:

    FrameList();
    ~FrameList();

    Frame *ActiveReference();
    void SetActiveRef(int);
    int AddRefFrame(Frame *, char *, AmberParm *,int,std::string&);
    int AddFrame(Frame *, AmberParm *);
    AmberParm *GetFrameParm(int);
    int GetFrameIndex(char *);
    int GetFrameIndexByTag(std::string &);
    Frame *GetFrame(int idx);
    int ReplaceFrame(int, Frame *, AmberParm *);
    void Info();
    const char *FrameName(int);

    int NumFrames() { return Nframe; }
};
#endif
