#ifndef INC_IMAGE_LIST_UNIT_FIRST_H
#define INC_IMAGE_LIST_UNIT_FIRST_H
#include "Image_List_Unit.h"
namespace Image {
class List_Unit_First : public List_Unit {
  public:
    List_Unit_First() {}
    Vec3 GetCoord(unsigned int idx, Frame const& frm) const { return frm.XYZ(units_[idx].Front()); }
};
}
#endif
