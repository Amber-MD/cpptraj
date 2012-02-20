// FrameList
#include "FrameList.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
FrameList::FrameList() {
  Nframe=0;
  referenceFrameNum = 0;
  //firstFrameNum = -1;
}

// DESTRUCTOR
FrameList::~FrameList() {
  for (int i=0; i<Nframe; i++)
    delete frameList[i];
}

// FrameList::ActiveReference()
/** Return the address of the frame pointed to by referenceFrameNum.
  */
Frame *FrameList::ActiveReference() {
  if (Nframe<1) return NULL;
  return frameList[referenceFrameNum];
}

// FrameList::SetActiveRef()
/** Set the given frame list number as the active reference.
  */
void FrameList::SetActiveRef(int numIn) {
  if (numIn < 0 || numIn >= Nframe) {
    mprintf("Warning: FrameList::SetActiveRef: Ref # %i out of bounds.\n",numIn);
    return;
  }
  referenceFrameNum = numIn;
}

// FrameList::AddRefFrame()
/** Add given Frame to the FrameList. Store trajectory name and frame number 
  * that this frame came from in frameNames and frameNums respectively. Store 
  * the associated parm in FrameParm. 
  */
int FrameList::AddRefFrame(Frame *F, char *name, AmberParm *P, int framenum,
                           std::string &RefTag) 
{
  std::string nameCopy;

  if (F==NULL || name==NULL) return 1;

  // Check that name and RefTag are not currently in use
  if (GetFrameIndex(name)!=-1) {
    mprintf("Warning: Reference with name %s already exists!\n",name);
    //return 1;
  }
  if (GetFrameIndexByTag(RefTag)!=-1) {
    mprintf("Warning: Reference with tag %s already exists!\n",RefTag.c_str());
    //return 1;
  }

  frameList.push_back(F);
  nameCopy.assign(name);
  frameNames.push_back(nameCopy);
  frameTags.push_back(RefTag);
  frameNums.push_back(framenum);
  FrameParm.AddParm(P);
  Nframe++;
  return 0;
}

// FrameList::AddFrame()
/** Add given Frame to the FrameList. Store the associated parm in FrameParm.
  */
int FrameList::AddFrame(Frame *F, AmberParm *P) {
  if (F==NULL || P==NULL) return 1;
  frameList.push_back(F);
  FrameParm.AddParm(P);
  Nframe++;
  return 0;
}

// FrameList::GetFrameIndex()
/** Return index of frame in the frame list specified by name.
  */
int FrameList::GetFrameIndex(char *name) {
  int idx;
  std::string RefTag;

  // if first char of name is a bracket, assume tag.
  if (name[0]=='[') {
    RefTag.assign( name );
    idx = GetFrameIndexByTag( RefTag );
  
  // Otherwise assume filename or base filename
  } else {
    idx = -1;
    int fn_end = (int) frameNames.size();
    for (int fn=0; fn < fn_end; fn++) {
      if ( frameNames[fn].compare(name)==0 ) { idx=fn; break; }
    }
  }
  
  return idx;
}

// FrameList::GetFrameIndexByTag()
/** Return the index of frame in frame list specified by tag.
  */
int FrameList::GetFrameIndexByTag(std::string &RefTag) {
  int fn;
  if (RefTag.empty()) return -1;
  int fn_end = (int) frameTags.size();
  for (fn=0; fn < fn_end; fn++) 
    if (frameTags[fn].compare( RefTag )==0) return fn;
  return -1;
}

// FrameList::GetFrameParm()
/** Given index of frame, return parm in FrameParm
  */
AmberParm *FrameList::GetFrameParm(int idx) {
  return FrameParm.GetParm(idx);
}

// FrameList::GetFrame()
/** Return the frame in the frame list specified by index.
  */
Frame *FrameList::GetFrame(int idx) {
  if (idx<0 || idx>=Nframe) return NULL;
  return frameList[idx];
}

// FrameList::ReplaceFrame()
/** Replace the frame/parm at the given position with the given frame/parm.
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

// FrameList::Info()
/** Print a list of trajectory names that frames have been taken from.
  */
void FrameList::Info() {
  if (Nframe==0) {
    mprintf("  No frames defined.\n");
    return;
  }
  if (!frameNames.empty()) {
    mprintf("  The following %i frames have been defined:\n",Nframe);
    for (int fn=0; fn < Nframe; fn++) { 
      if (!frameTags[fn].empty())
        mprintf("    %i: %s frame %i\n",fn,frameTags[fn].c_str(),
                frameNums[fn]+OUTPUTFRAMESHIFT);
      else
        mprintf("    %i: %s frame %i\n",fn,frameNames[fn].c_str(),
                frameNums[fn]+OUTPUTFRAMESHIFT);
    }
  } else {
    mprintf("  %i frames have been defined.\n",Nframe);
  }
  mprintf("\tActive reference frame for masks is %i\n",referenceFrameNum);
}

// FrameList::FrameName()
/** Return name of given frame.
  */
const char *FrameList::FrameName(int idx) {
  if (idx<0 || idx>=(int)frameNames.size()) return NULL;
  return frameNames[idx].c_str();
}

