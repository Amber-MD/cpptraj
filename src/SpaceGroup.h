#ifndef INC_SPACEGROUP_H
#define INC_SPACEGROUP_H
#include <string>
#include <vector>
#include <map>
#include "Matrix_3x3.h"
#include "Vec3.h"
class SpaceGroup {
  public:
    SpaceGroup();
    /// \return a number indicating the space group that corresponds to given text.
    int ID(std::string const&);
    /// Create a list of symmetry operations based on current space group.
    int LoadSymmOps(int, int, int, std::vector<Matrix_3x3>&, std::vector<Vec3>&) const;
  private:
    typedef std::map<std::string, int> MapType;
    typedef std::pair<std::string, int> PairType;
    MapType IdMap_; ///< Map space group strings to internal IDs
    int sgID_; ///< Current space group ID
};
#endif
