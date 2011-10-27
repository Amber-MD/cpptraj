// FrameList
#include "FrameList.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
FrameList::FrameList() {
  Nframe=0;
}

// DESTRUCTOR
FrameList::~FrameList() {
  for (int i=0; i<Nframe; i++)
    delete frameList[i];
}

/* FrameList::AddRefFrame()
 * Add given Frame to the FrameList. Store trajectory name and frame number 
 * that this frame came from in frameNames and frameNums respectively. Store 
 * the associated parm in FrameParm. 
 */
int FrameList::AddRefFrame(Frame *F, char *name, AmberParm *P, int framenum) {
  std::string nameCopy;

  if (F==NULL || name==NULL) return 1;

  frameList.push_back(F);
  nameCopy.assign(name);
  frameNames.push_back(nameCopy);
  frameNums.push_back(framenum);
  FrameParm.Add(P);
  Nframe++;
  return 0;
}

/* FrameList::AddFrame()
 * Add given Frame to the FrameList. Store the associated parm in FrameParm.
 */
int FrameList::AddFrame(Frame *F, AmberParm *P) {
  if (F==NULL || P==NULL) return 1;
  frameList.push_back(F);
  FrameParm.Add(P);
  Nframe++;
  return 0;
}

/* FrameList::GetFrameIndex()
 * Return index of frame in the frame list specified by name.
 */
int FrameList::GetFrameIndex(char *name) {
  int idx = -1;
  int fn_end = (int) frameNames.size();
  for (int fn=0; fn < fn_end; fn++) {
    if ( frameNames[fn].compare(name)==0 ) { idx=fn; break; }
  }
  
  return idx;
}

/* FrameList::GetFrameParm()
 * Given index of frame, return parm in FrameParm
 */
AmberParm *FrameList::GetFrameParm(int idx) {
  return FrameParm.GetParm(idx);
}

/* FrameList::GetFrame()
 * Return the frame in the frame list specified by index.
 */
Frame *FrameList::GetFrame(int idx) {
  if (idx<0 || idx>=Nframe) return NULL;
  return frameList[idx];
}

/* FrameList::Replace()
 * Replace the frame/parm at the given position with the given frame/parm.
 * The old frame is deleted. 
 */
int FrameList::Replace(int idx, Frame *newFrame, AmberParm *newParm) {
  if (newFrame==NULL || newParm==NULL) return 1;
  if (idx<0 || idx>=Nframe) return 1;
  delete frameList[idx];
  frameList[idx]=newFrame;
  if (FrameParm.Replace(idx,newParm)) return 1;
  return 0;
}

/* FrameList::Info()
 * Print a list of trajectory names that frames have been taken from.
 */
void FrameList::Info() {
  if (Nframe==0) {
    mprintf("  No frames defined.\n");
    return;
  }
  if (!frameNames.empty()) {
    mprintf("  The following %i frames have been defined:\n",Nframe);
    for (int fn=0; fn < Nframe; fn++) 
      mprintf("    %i: %s frame %i\n",fn,frameNames[fn].c_str(),
              frameNums[fn]+OUTPUTFRAMESHIFT);
  } else {
    mprintf("  %i frames have been defined.\n",Nframe);
  }
}

/* FrameList::FrameName()
 * Return name of given frame.
 */
const char *FrameList::FrameName(int idx) {
  if (idx<0 || idx>=(int)frameNames.size()) return NULL;
  return frameNames[idx].c_str();
}

