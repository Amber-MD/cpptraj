#ifndef INC_IMAGE_LIST_UNIT_H
#define INC_IMAGE_LIST_UNIT_H
#include "Image_List.h"
#include "Frame.h"
#include <vector>
class Unit;
namespace Image {
/// List with entities defined by Units
class List_Unit : public List {
  public:
    List_Unit() {}
    virtual ~List_Unit() {}
    unsigned int nEntities() const { return units_.size(); }
    int SetupList(Topology const&, std::string const&);
    void DoTranslation(Frame& frm, unsigned int idx, Vec3 const& boxTrans) const {
      frm.Translate(boxTrans, units_[idx]);
    }
    void CopyFrame(Frame& refIn, unsigned int idx, Frame const& tgtIn) const {
      refIn.CopyFrom(tgtIn, units_[idx]);
    }

    void AddUnit(Unit const&);
  protected:
    typedef std::vector<Unit> Uarray;
    Uarray units_;
};
}
#endif
