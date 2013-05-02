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
    static void CheckEven(int&, const char*);
    static const char* HelpText;
    static DataSet_GridFlt* AllocateGrid(DataSetList&,std::string const&,
                                         int,int,int,
                                         double,double,double,
                                         double,double,double);
    DataSet_GridFlt* GridInit(const char*, ArgList&, DataSetList&);
    void GridInfo(DataSet_GridFlt const&);
    int GridSetup(Topology const&);
    void GridFrame(Frame const&, AtomMask const&, DataSet_GridFlt&);
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
  Vec3 center;
  if      (mode_==BOX)    center = currentFrame.BoxCrd().Center(); 
  else if (mode_==CENTER) center = currentFrame.VGeometricCenter( centerMask_ );
  else                    center.Zero(); // mode_==ORIGIN
  for (AtomMask::const_iterator atom = mask.begin(); atom != mask.end(); ++atom)
    grid.Increment( Vec3(currentFrame.XYZ(*atom)) - center, increment_ );
}
#endif
