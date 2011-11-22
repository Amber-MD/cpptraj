// FrameList
#include "FrameList.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
FrameList::FrameList() {
  Nframe=0;
  referenceFrameNum = 0;
  firstFrameNum = -1;
}

// DESTRUCTOR
FrameList::~FrameList() {
  for (int i=0; i<Nframe; i++)
    delete frameList[i];
}

/* FrameList::ActiveRefCoords()
 * Return the coords of frame pointed to by referenceFrameNum.
 */
double *FrameList::ActiveReference() {
  if (Nframe<1) return NULL;
  return frameList[referenceFrameNum]->X;
}

/* FrameList::SetActiveRef()
 * Set the given frame list number as the active reference.
 */
void FrameList::SetActiveRef(int numIn) {
  if (numIn < 0 || numIn >= Nframe) {
    mprintf("Warning: FrameList::SetActiveRef: Ref # %i out of bounds.\n",numIn);
    return;
  }
  referenceFrameNum = numIn;
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
  FrameParm.AddParm(P);
  Nframe++;
  return 0;
}

// FrameList::AddFirstFrame()
/// Add the input frame as the first frame of the trajectory.
/** This is a special case of AddRefFrame. Actions like Rmsd may
  * require the first frame specifically. This will be called by
  * the first such action, subsequent actions will get this frame.
  * The frame will be given the name of the parm.
  */
int FrameList::AddFirstFrame(Frame *frameIn, AmberParm *parmIn) {
  if (AddRefFrame(frameIn, parmIn->parmName, parmIn, 0)) return 1;
  // Since Nframe has now been incremented, actual position in 
  // frameList is Nframe - 1.
  firstFrameNum = Nframe - 1;
  return 0;
}

// FrameList::GetFirstFrame()
/// Return the frame indicated as the first frame if previously set up.
/** If a frame has been added with the AddFirstFrame command return it,
  * otherwise return NULL.
  */
Frame *FrameList::GetFirstFrame() {
  if (firstFrameNum!=-1) return frameList[firstFrameNum];
  return NULL;
}

// FrameList::GetFirstFrameParm()
AmberParm *FrameList::GetFirstFrameParm() {
  if (firstFrameNum!=-1) return FrameParm.GetParm(firstFrameNum);
  return NULL;
}

/* FrameList::AddFrame()
 * Add given Frame to the FrameList. Store the associated parm in FrameParm.
 */
int FrameList::AddFrame(Frame *F, AmberParm *P) {
  if (F==NULL || P==NULL) return 1;
  frameList.push_back(F);
  FrameParm.AddParm(P);
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

/* FrameList::ReplaceFrame()
 * Replace the frame/parm at the given position with the given frame/parm.
 * The old frame is deleted. 
 */
int FrameList::ReplaceFrame(int idx, Frame *newFrame, AmberParm *newParm) {
  if (newFrame==NULL || newParm==NULL) return 1;
  if (idx<0 || idx>=Nframe) return 1;
  delete frameList[idx];
  frameList[idx]=newFrame;
  if (FrameParm.ReplaceParm(idx,newParm)) return 1;
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
  mprintf("\tActive reference frame for masks is %i\n",referenceFrameNum);
}

/* FrameList::FrameName()
 * Return name of given frame.
 */
const char *FrameList::FrameName(int idx) {
  if (idx<0 || idx>=(int)frameNames.size()) return NULL;
  return frameNames[idx].c_str();
}

