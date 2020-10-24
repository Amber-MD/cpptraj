#ifndef INC_IMAGE_LIST_MASK_H
#define INC_IMAGE_LIST_MASK_H
#include "Image_List.h"
#include "Frame.h"
#include "AtomMask.h"
namespace Image {
/// List with entities defined by atom mask 
class List_Mask : public List {
  public:
    List_Mask() {}
    virtual ~List_Mask() {}
    unsigned int nEntities() const { return mask_.Nselected(); }
    int SetupList(Topology const&, std::string const&);
    Vec3 GetCoord(unsigned int idx, Frame const& frm) const { return frm.XYZ(mask_[idx]); }
    void DoTranslation(Frame& frm, unsigned int idx, Vec3 const& boxTrans) const {
      frm.Translate(boxTrans, mask_[idx]);
    }
    Unit AllEntities() const;
    void PrintEntities() const;
  protected:
    AtomMask mask_;
};
}
#endif
