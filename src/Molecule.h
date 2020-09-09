#ifndef INC_MOLECULE_H
#define INC_MOLECULE_H
#include "Unit.h"
// Class: Molecule
/// Hold information for a molecule
class Molecule {
  public:
    Molecule() : isSolvent_(false) {}
    /// CONSTRUCTOR - Single segment from beg to end
    Molecule(int beg, int end) : isSolvent_(false), unit_(Unit(beg, end)) {}

    void SetSolvent()   { isSolvent_ = true;  }
    void SetNoSolvent() { isSolvent_ = false; }
    Unit& ModifyUnit()  { return unit_;       }

    bool IsSolvent()        const { return isSolvent_;   }
    unsigned int NumAtoms() const { return unit_.Size(); }
    Unit const& MolUnit()   const { return unit_; }

  private:
    bool isSolvent_;
    Unit unit_;      ///< Hold all molecule segments.
};
#endif
