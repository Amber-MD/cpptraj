#ifndef INC_GRIDACTION_H
#define INC_GRIDACTION_H
#include "ArgList.h"
#include "DataSetList.h"
#include "DataSet_GridFlt.h"
#include "Topology.h"
/// Class for setting up a grid within an action.
class GridAction {
  public:
    /// Indicate which kind of gridding to perform
    enum GridModeType { ORIGIN = 0, BOX, CENTER };
    GridAction() {}
    static const char* HelpText;
    DataSet_GridFlt* GridInit(const char*, ArgList&, DataSetList&);
    void GridInfo(DataSet_GridFlt const&);
    int GridSetup(Topology const&);
    GridModeType GridMode()      const { return mode_;       }
    AtomMask const& CenterMask() const { return centerMask_; }
    float Increment()            const { return increment_;  }
  private:
    GridModeType mode_;
    AtomMask centerMask_;
    float increment_;     ///< Set to -1 if negative, 1 if not.
};
#endif
