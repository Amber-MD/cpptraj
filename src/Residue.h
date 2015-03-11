#ifndef INC_RESIDUE_H
#define INC_RESIDUE_H
#include "NameType.h"
// Class: Residue
/// Hold information for a residue.
class Residue {
  public:
    Residue() : resname_(""), firstAtom_(0), lastAtom_(0), icode_(' ') {}
    Residue(int onum, NameType const& resname, int firstAtomIn) :
      resname_(resname), firstAtom_(firstAtomIn), originalResNum_(onum), icode_(' ')
    {}
    inline void SetLastAtom(int i)      { lastAtom_ = i;          }
    inline void SetOriginalNum(int i)   { originalResNum_ = i;    }
    inline void SetIcode(char c)        { icode_ = c;             }
    /// \return First atom in residue, indexing from 0
    inline int FirstAtom()        const { return firstAtom_;      }
    /// \return Atom _after_ the last in residue, indexing from 0
    inline int LastAtom()         const { return lastAtom_;       }
    inline int OriginalResNum()   const { return originalResNum_; }
    inline char Icode()           const { return icode_;          }
    inline const char *c_str()    const { return *resname_;       }
    inline NameType const& Name() const { return resname_;        }
    inline int NumAtoms()         const { return (lastAtom_ - firstAtom_); }
    inline bool NameIsSolvent()   const {
      return (resname_=="WAT " || resname_=="HOH " || resname_=="TIP3");
    }
    /// Convert 3-letter residue code to single letter.
    static char ConvertResName(std::string const&);
    /// Convert this residue name to single letter.
    char SingleCharName() const { return ConvertResName( *resname_ ); }
  private:
    NameType resname_;   ///< Residue name.
    int firstAtom_;      ///< Index of first atom (from 0).
    int lastAtom_;       ///< Atom index after last atom in residue.
    int originalResNum_; ///< Original residue number.
    char icode_;         ///< Residue insertion code.
};
#endif
