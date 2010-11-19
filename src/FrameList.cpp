// FrameList
#include <cstdio>
#include <cstdlib>
#include "FrameList.h"

// CONSTRUCTOR
FrameList::FrameList() {
  frameList=NULL;
  Nframe=0;
}

// DESTRUCTOR
FrameList::~FrameList() {
  int i;

  if (frameList!=NULL) {
    for (i=0; i<Nframe; i++)
      delete frameList[i];
    free(frameList);
  }
}

/* FrameList::Add()
 * Add given Frame to the FrameList. Store trajectory name that this frame
 * came from in frameNames. Store the associated parm in FrameParm. 
 */
int FrameList::Add(Frame *F, char *name, AmberParm *P) {

  if (F==NULL || name==NULL) return 1;

  frameList=(Frame**) realloc(frameList, (Nframe+1) * sizeof(Frame*));
  frameList[Nframe]=F;
  frameNames.Add(name);
  Nframe++;
  FrameParm.Add(P);
  return 0;
}

/* FrameList::GetFrameIndex()
 * Return index of frame in the frame list specified by name.
 */
//Frame *FrameList::GetFrameIndex(char *name) {
int FrameList::GetFrameIndex(char *name) {
  //int idx;
  
  return frameNames.getKeyIndex(name);
  /*if (idx==-1) {
    fprintf(stdout,"Could not find frame %s in frame list.\n",name);
    return NULL;
  } else
    return frameList[idx];*/
}

/*
 * FrameList::GetFrameParm()
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

/* FrameList::Info()
 * Print a list of trajectory names that frames have been taken from.
 */
void FrameList::Info() {
  if (Nframe==0) {
    fprintf(stdout,"  No frames defined.\n");
    return;
  }
  fprintf(stdout,"  The following %i frames have been defined:\n",Nframe);
  frameNames.print();
}
