#ifndef INC_IMAGE_LIST_UNIT_COM_H
#define INC_IMAGE_LIST_UNIT_COM_H
#include "Image_List_Unit.h"
namespace Image {
class List_Unit_CoM : public List_Unit {
  public:
    List_Unit_CoM() {}
    Vec3 GetCoord(unsigned int idx, Frame const& frm) const { return frm.VCenterOfMass(units_[idx]); }
};
}
#endif
