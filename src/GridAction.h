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
    enum GridModeType { ORIGIN = 0, BOX, MASKCENTER, SPECIFIEDCENTER };
    GridAction() : increment_(1.0) {}
    static const char* HelpText;
    DataSet_GridFlt* GridInit(const char*, ArgList&, DataSetList&);
    void GridInfo(DataSet_GridFlt const&);
    int GridSetup(Topology const&);
    inline void GridFrame(Frame const&, AtomMask const&, DataSet_GridFlt&);
    GridModeType GridMode()      const { return mode_;       }
    AtomMask const& CenterMask() const { return centerMask_; }
    float Increment()            const { return increment_;  }
  private:
    GridModeType mode_;
    AtomMask centerMask_;
    float increment_;     ///< Set to -1 if negative, 1 if not.
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
void GridAction::GridFrame(Frame const& currentFrame, AtomMask const& mask, 
                           DataSet_GridFlt& grid) 
{
  Vec3 offset;
  if      (mode_==BOX)        offset = currentFrame.BoxCrd().Center(); 
  else if (mode_==MASKCENTER) offset = currentFrame.VGeometricCenter( centerMask_ );
  else                        offset.Zero(); // mode_==ORIGIN/SPECIFIEDCENTER
  for (AtomMask::const_iterator atom = mask.begin(); atom != mask.end(); ++atom)
    grid.Increment( Vec3(currentFrame.XYZ(*atom)) - offset, increment_ );
}
#endif
