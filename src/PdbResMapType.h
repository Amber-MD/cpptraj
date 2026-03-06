#ifndef INC_PDBRESMAPTYPE_H
#define INC_PDBRESMAPTYPE_H
#include <string>
#include "NameType.h"
#include "Structure/StructureEnum.h"
namespace Cpptraj {
/// Hold PDB residue name, corresponding unit name and residue terminal type
class PdbResMapType {
  public:
    //PdbResMapType() : termType_(Structure::NON_TERMINAL) {}
    PdbResMapType(std::string const& unitIn, NameType const& nameIn, Structure::TerminalType termIn) :
      unitName_(unitIn), pdbName_(nameIn), termType_(termIn) {}

    std::string const& UnitName() const { return unitName_; }
    NameType const& PdbName() const { return pdbName_; }
    Structure::TerminalType TermType() const { return termType_; }
  private:
    std::string unitName_;
    NameType pdbName_;
    Structure::TerminalType termType_;
};
}
#endif
