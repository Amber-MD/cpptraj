#ifndef INC_FRAMELIST_H
#define INC_FRAMELIST_H
#include "ArgList.h"
#include "Frame.h"
#include "Topology.h"
#include "FileList.h"
// Class: FrameList
/// Hold a series of Frame classes. 
/** Optionally hold a corresponding name and number (in the case of 
  * reference frames, this is the name of the reference coordinate 
  * file and the frame number). 
  */
// TODO: Eventually store a vector of Frames, not Frame*s
class FrameList : public FileList {
  public:
    FrameList();
    ~FrameList();
    static void Help();
    Frame *ActiveReference();
    void SetActiveRef(int);
    int AddReference(ArgList *, Topology *);
    //int AddRefFrame(Frame *, char *, const char *,Topology *,int,std::string&);
    int AddFrame(Frame *, Topology *);
    Topology *GetFrameParm(int);
    Frame *GetFrame(int idx);
    int ReplaceFrame(int, Frame *, Topology *);
    void Info();
    const char *FrameName(int);

    inline int NumFrames() { 
      return (int)frames_.size();
    }
  private:
    std::vector<Frame*> frames_;
    std::vector<Topology*> parms_;
    std::vector<int> nums_;
    std::vector<Topology*> StrippedRefParms_;
    int refFrameNum_;
};
#endif
