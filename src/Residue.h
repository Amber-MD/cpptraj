#ifndef INC_RESIDUE_H
#define INC_RESIDUE_H
#include "NameType.h"
// Class: Residue
/// Hold information for a residue.
class Residue {
  public:
    Residue() : firstAtom_(0), lastAtom_(0), resname_("") {}
    Residue(int onum, NameType const& resname, int firstAtomIn) :
      firstAtom_(firstAtomIn), originalResNum_(onum), resname_(resname)
    {}
    Residue(int onum, NameType const& resname, int firstAtomIn, int lastAtomIn) :
      firstAtom_(firstAtomIn), lastAtom_(lastAtomIn), originalResNum_(onum), 
      resname_(resname)
    {}
    inline void SetLastAtom(int i)      { lastAtom_ = i;          }
    inline void SetOriginalNum(int i)   { originalResNum_ = i;    }
    /// \return First atom in residue, indexing from 0
    inline int FirstAtom()        const { return firstAtom_;      }
    /// \return Atom _after_ the last in residue, indexing from 0
    inline int LastAtom()         const { return lastAtom_;       }
    inline int OriginalResNum()   const { return originalResNum_; }
    inline const char *c_str()    const { return *resname_;       }
    inline NameType const& Name() const { return resname_;        }
    inline int NumAtoms()         const { return (lastAtom_ - firstAtom_); }
  private:
    /// The first atom in the residue, atom numbering starts from 0
    int firstAtom_;
    /// Actually the atom _after_ the last atom in the residue.
    int lastAtom_;
    /// The original residue number.
    int originalResNum_;
    NameType resname_;
};
#endif
