#ifndef INC_RESIDUE_H
#define INC_RESIDUE_H
#include "NameType.h"
// Class: Residue
/// Hold information for a residue.
class Residue {
  public:
    Residue() : firstAtom_(0), lastAtom_(0), resname_("") {}
    Residue(NameType const& resname, int firstAtomIn) :
      firstAtom_(firstAtomIn), resname_(resname)
    {}
    Residue(NameType const& resname, int firstAtomIn, int lastAtomIn) :
      firstAtom_(firstAtomIn), lastAtom_(lastAtomIn), resname_(resname)
    {}
    inline void SetLastAtom(int i)      { lastAtom_ = i;     }
    inline int FirstAtom()        const { return firstAtom_; }
    inline int LastAtom()         const { return lastAtom_;  }
    inline const char *c_str()    const { return *resname_;  }
    inline NameType const& Name() const { return resname_;   }
    inline int NumAtoms()         const { return (lastAtom_ - firstAtom_); }
  private:
    int firstAtom_;
    int lastAtom_;
    NameType resname_;
};
#endif
