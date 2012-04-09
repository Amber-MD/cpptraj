// FrameList
#include "FrameList.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
FrameList::FrameList() : 
  refFrameNum_(0)
{}

// DESTRUCTOR
FrameList::~FrameList() {
  for (std::vector<Frame*>::iterator frame = frames_.begin();
                                     frame != frames_.end(); frame++)
    delete *frame;
}

// FrameList::ActiveReference()
/** Return the address of the frame pointed to by refFrameNum_.
  */
Frame *FrameList::ActiveReference() {
  if (frames_.empty()) return NULL;
  return frames_[refFrameNum_];
}

// FrameList::SetActiveRef()
/** Set the given frame list number as the active reference.
  */
void FrameList::SetActiveRef(int numIn) {
  if (numIn < 0 || numIn >= (int)frames_.size()) {
    mprintf("Warning: FrameList::SetActiveRef: Ref # %i out of bounds.\n",numIn);
    return;
  }
  refFrameNum_ = numIn;
}

// FrameList::AddRefFrame()
/** Add given Frame to the FrameList. Store trajectory name and frame number 
  * that this frame came from in frameNames and frameNums respectively. Store 
  * the associated parm in FrameParm. 
  */
int FrameList::AddRefFrame(Frame *F, char *name, const char *filename, Topology *P, 
                           int framenum, std::string &RefTag) 
{
  if (F==NULL || name==NULL) return 1;

  // Check that name and RefTag are not currently in use
  if (FindName(name)!=-1) {
    mprintf("Warning: Reference with name %s already exists!\n",name);
    //return 1;
  }
  if (FindName(RefTag)!=-1) {
    mprintf("Warning: Reference with tag %s already exists!\n",RefTag.c_str());
    //return 1;
  }

  frames_.push_back(F);
  AddNames(name, filename, RefTag);
  nums_.push_back(framenum);
  parms_.push_back(P);
  return 0;
}

// FrameList::AddFrame()
/** Add given Frame to the FrameList. Store the associated parm in FrameParm.
  */
int FrameList::AddFrame(Frame *F, Topology *P) {
  if (F==NULL || P==NULL) return 1;
  frames_.push_back(F);
  parms_.push_back(P);
  return 0;
}

// FrameList::GetFrameParm()
/** Given index of frame, return parm in FrameParm
  */
Topology *FrameList::GetFrameParm(int idx) {
  return parms_[idx];
}

// FrameList::GetFrame()
/** Return the frame in the frame list specified by index.
  */
Frame *FrameList::GetFrame(int idx) {
  if (idx<0 || idx>=(int)frames_.size()) return NULL;
  return frames_[idx];
}

// FrameList::ReplaceFrame()
/** Replace the frame/parm at the given position with the given frame/parm.
  * The old frame is deleted. 
  */
int FrameList::ReplaceFrame(int idx, Frame *newFrame, Topology *newParm) {
  if (newFrame==NULL || newParm==NULL) return 1;
  if (idx<0 || idx>=(int)frames_.size()) return 1;
  delete frames_[idx];
  frames_[idx] = newFrame;
  parms_[idx] = newParm;
  return 0;
}

// FrameList::Info()
/** Print a list of trajectory names that frames have been taken from.
  */
void FrameList::Info() {
  if (frames_.empty()) {
    mprintf("  No frames defined.\n");
    return;
  }
  if (HasNames()) {
    mprintf("  The following %zu frames have been defined:\n",frames_.size());
    for (int fn=0; fn < (int)frames_.size(); fn++) { 
      if (!Tag(fn).empty())
        mprintf("    %i: %s frame %i\n",fn,Tag(fn).c_str(),
                nums_[fn]+OUTPUTFRAMESHIFT);
      else
        mprintf("    %i: %s frame %i\n",fn,Name(fn).c_str(),
                nums_[fn]+OUTPUTFRAMESHIFT);
    }
  } else {
    mprintf("  %zu frames have been defined.\n",frames_.size());
  }
  mprintf("\tActive reference frame for masks is %i\n",refFrameNum_);
}

// FrameList::FrameName()
/** Return name of given frame.
  */
const char *FrameList::FrameName(int idx) {
  if (idx<0 || idx>=(int)frames_.size()) return NULL;
  if (Name(idx).empty()) return NULL;
  return Name(idx).c_str();
}

