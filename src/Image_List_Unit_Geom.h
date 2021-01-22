#ifndef INC_IMAGE_LIST_UNIT_GEOM_H
#define INC_IMAGE_LIST_UNIT_GEOM_H
#include "Image_List_Unit.h"
namespace Image {
class List_Unit_Geom : public List_Unit {
  public:
    List_Unit_Geom() {}
    Vec3 GetCoord(unsigned int idx, Frame const& frm) const { return frm.VGeometricCenter(units_[idx]); }
};
}
#endif
