#ifndef INC_TRAJINLIST_H
#define INC_TRAJINLIST_H
/// Class: TrajinList
#include "CoordFileList.h"
class TrajinList : public CoordFileList {
  public:
    TrajinList();
    ~TrajinList();
    // Inherited: Add a traj file to the list based on input from arg list
    int AddTrajin(char*, ArgList *, AmberParm *);
    // TRAJIN: Set up frames to be processed 
    int SetupFrames();
};
#endif

