#ifndef INC_RESIDUE_H
#define INC_RESIDUE_H
#include "NameType.h"
// Class: Residue
/// Hold information for a residue.
class Residue {
  public:
    Residue() : firstAtom_(0), resname_("") {}
    Residue(NameType const& resname, int firstAtomIn) :
      firstAtom_(firstAtomIn),
      resname_(resname)
    {}

    inline int FirstAtom() const        { return firstAtom_; }
    inline const char *c_str() const    { return *resname_;  }
    inline NameType const& Name() const { return resname_;   }
  private:
    int firstAtom_;
    NameType resname_;
};
#endif
