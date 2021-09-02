#ifndef INC_IMAGE_LIST_PAIR_COM_H
#define INC_IMAGE_LIST_PAIR_COM_H
#include "Image_List_Pair.h"
namespace Image {
class List_Pair_CoM : public List_Pair {
  public:
    List_Pair_CoM() {}
    Vec3 GetCoord(unsigned int idx, Frame const& frm) const { return frm.VCenterOfMass(begin_[idx], end_[idx]); }
};
}
#endif
