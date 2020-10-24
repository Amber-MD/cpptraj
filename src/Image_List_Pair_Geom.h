#ifndef INC_IMAGE_LIST_PAIR_GEOM_H
#define INC_IMAGE_LIST_PAIR_GEOM_H
#include "Image_List_Pair.h"
namespace Image {
class List_Pair_Geom : public List_Pair {
  public:
    List_Pair_Geom() {}
    Vec3 GetCoord(unsigned int idx, Frame const& frm) const { return frm.VGeometricCenter(begin_[idx], end_[idx]); }
};
}
#endif
