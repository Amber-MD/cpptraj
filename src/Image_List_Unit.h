#ifndef INC_IMAGE_LIST_UNIT_H
#define INC_IMAGE_LIST_UNIT_H
#include "Image_List.h"
#include "Frame.h"
#include <vector>
class Unit;
namespace Image {
/// List with entities defined by Units
class List_Unit : public List {
    typedef std::vector<Unit> Uarray;
  public:
    List_Unit() {}
    virtual ~List_Unit() {}
    unsigned int nEntities() const { return units_.size(); }
    int SetupList(Topology const&, std::string const&);
    void DoTranslation(Frame& frm, unsigned int idx, Vec3 const& boxTrans) const {
      frm.Translate(boxTrans, units_[idx]);
    }
    Unit AllEntities() const;
    void PrintEntities() const;

    void AddUnit(Unit const&);
    typedef Uarray::const_iterator const_iterator;
    const_iterator begin() const { return units_.begin(); }
    const_iterator end()   const { return units_.end(); }
    bool empty()           const { return units_.empty(); }
  protected:
    Uarray units_;
};
}
#endif
