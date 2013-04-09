#ifndef INC_GRIDACTION_H
#define INC_GRIDACTION_H
#include "ArgList.h"
#include "DataSetList.h"
#include "DataSet_3D.h"
/// Class for setting up a grid within an action.
class GridAction {
  public:
    GridAction() {}
    static const char* HelpText;
    DataSet_3D* GridInit(const char*, ArgList&, DataSetList&);
};
#endif
