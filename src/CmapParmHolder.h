#ifndef INC_PARM_CMAPPARMHOLDER_H
#define INC_PARM_CMAPPARMHOLDER_H
class CmapGridType;
#include <vector>
#include "Parm/ParmEnum.h"
//namespace Cpptraj {
//namespace Parm {
/// Hold CMAP terms
class CmapParmHolder {
    typedef std::vector<CmapGridType> Carray;
  public:
    CmapParmHolder();

    /// const iterator
    typedef Carray::const_iterator const_iterator;
    /// const begin iterator
    const_iterator begin() const { return CMAP_.begin(); }
    /// const end iterator
    const_iterator end() const { return CMAP_.end(); }
    /// \return Number of CMAP terms
    Carray::size_type size() const { return CMAP_.size(); }
    /// \return true if no CMAP terms
    bool empty() const { return CMAP_.empty(); }
    /// \return CMAP at index
    CmapGridType const& operator[](int idx) const { return CMAP_[idx]; }

    /// \return Modifiable CMAP at index
    
    /// Add/update cmap term; optionally allow update, pass in debug level
    Cpptraj::Parm::RetType AddParm(CmapGridType const&, bool, int);
  private:
    Carray CMAP_;
};
//}
//}
#endif
