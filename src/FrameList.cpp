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

/* FrameList::Add()
 * Add given Frame to the FrameList. Store trajectory name that this frame
 * came from in frameNames. Store the associated parm in FrameParm. 
 */
int FrameList::Add(Frame *F, char *name, AmberParm *P, int framenum) {
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

/* 
 * FrameList::GetFrameIndex()
 * Return index of frame in the frame list specified by name.
 */
int FrameList::GetFrameIndex(char *name) {
  int idx = -1;

  for (int fn=0; fn < Nframe; fn++) {
    if ( frameNames[fn].compare(name)==0 ) { idx=fn; break; }
  }
  
  return idx;
}

/*
 * FrameList::GetFrameParm()
 * Given index of frame, return parm in FrameParm
 */
AmberParm *FrameList::GetFrameParm(int idx) {
  return FrameParm.GetParm(idx);
}

/* 
 * FrameList::GetFrame()
 * Return the frame in the frame list specified by index.
 */
Frame *FrameList::GetFrame(int idx) {
  if (idx<0 || idx>=Nframe) return NULL;
  return frameList[idx];
}

/* 
 * FrameList::Info()
 * Print a list of trajectory names that frames have been taken from.
 */
void FrameList::Info() {
  if (Nframe==0) {
    mprintf("  No frames defined.\n");
    return;
  }
  mprintf("  The following %i frames have been defined:\n",Nframe);
  for (int fn=0; fn < Nframe; fn++) 
    mprintf("    %i: %s frame %i\n",fn,frameNames[fn].c_str(),
            frameNums[fn]+OUTPUTFRAMESHIFT);
}
