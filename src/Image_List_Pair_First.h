#ifndef INC_IMAGE_LIST_PAIR_FIRST_H
#define INC_IMAGE_LIST_PAIR_FIRST_H
#include "Image_List_Pair.h"
namespace Image {
class List_Pair_First : public List_Pair {
  public:
    List_Pair_First() {}
    Vec3 GetCoord(unsigned int idx, Frame const& frm) const { return frm.XYZ(begin_[idx]); }
};
}
#endif
