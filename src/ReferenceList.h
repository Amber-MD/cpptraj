#ifndef INC_REFERENCELIST_H
#define INC_REFERENCELIST_H
#include "CoordFileList.h"
#include "FrameList.h"
// Class: ReferenceList
class ReferenceList : public CoordFileList {
    std::vector<bool> Average;
    std::vector<std::string> MaskExpressions;
    std::vector<std::string> RefTags;
    std::vector<AmberParm*> StrippedRefParms;
  public:
    
    ReferenceList();
    ~ReferenceList();
    // Add a traj file to the list with given access and associate with a parm
    int AddReference(char*,ArgList *A, AmberParm *);
    // REFERENCE: Set up frames to be processed
    int SetupRefFrames(FrameList *);
};
#endif

