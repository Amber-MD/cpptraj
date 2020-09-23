#ifndef INC_IMAGE_LIST_PAIR_H
#define INC_IMAGE_LIST_PAIR_H
#include "Image_List.h"
#include "Frame.h"
#include <vector>
namespace Image {
/// List with entities defined by a single range
class List_Pair : public List {
  public:
    List_Pair() {}
    virtual ~List_Pair() {}
    unsigned int nEntities() const { return begin_.size(); }
    int SetupList(Topology const&, std::string const&);
    void DoTranslation(Frame& frm, unsigned int idx, Vec3 const& boxTrans) const {
      frm.Translate(boxTrans, begin_[idx], end_[idx]);
    }
    Unit AllEntities() const;
    void PrintEntities() const;
  protected:
    typedef std::vector<int> Iarray;
    Iarray begin_;
    Iarray end_;
};
}
#endif
