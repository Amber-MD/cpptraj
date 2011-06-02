// TrajoutList
#ifndef INC_TRAJOUTLIST_H
#define INC_TRAJOUTLIST_H
#include "CoordFileList.h"

class TrajoutList : public CoordFileList {

  public:
    
    TrajoutList();
    ~TrajoutList();
    // Add a traj file to the list with given access and associate with a parm
    int Add(ArgList *A, AmberParm *);
    int Write(int , Frame *, AmberParm*);
    void Close();
};
#endif

