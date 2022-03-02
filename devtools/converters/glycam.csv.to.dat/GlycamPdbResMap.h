#ifndef INC_GLYCAMPDBRESMAP_H
#define INC_GLYCAMPDBRESMAP_H
#include <string>
#include <map>
#include "NameType.h"
class GlycamPdbResMap {
    typedef std::pair<NameType, char> PairType;
    typedef std::map<NameType, char> MapType;
  public:
    GlycamPdbResMap() {}
    int Load(std::string const&);
    unsigned int size() const { return pdb_to_glycam_.size(); }

    typedef MapType::const_iterator const_iterator;
    const_iterator begin() const { return pdb_to_glycam_.begin(); }
    const_iterator end()   const { return pdb_to_glycam_.end();   }
  private:
    MapType pdb_to_glycam_; ///< Map PDB residue names to Glycam 1 char names
};
#endif
