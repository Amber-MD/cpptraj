#ifndef INC_RESIDUE_H
#define INC_RESIDUE_H
#include "NameType.h"
// Class: Residue
/// Hold information for a residue.
class Residue {
  public:
    Residue() : resname_(""), firstAtom_(0), lastAtom_(0), icode_(' '), chainID_(' ') {}
    /// CONSTRUCTOR - Copy given Residue, set first and last atom indices.
    Residue(Residue const& r, int first, int last) :
      resname_(r.resname_), firstAtom_(first), lastAtom_(last),
      originalResNum_(r.originalResNum_), icode_(r.icode_), chainID_(r.chainID_)
    {}
    /// CONSTRUCTOR - Res name, original resnum, icode, chain ID
    Residue(NameType const& n, int r, char ic, char cid) :
            resname_(n), originalResNum_(r), icode_(ic), chainID_(cid) {}
    inline void SetFirstAtom(int i)        { firstAtom_ = i;      }
    inline void SetLastAtom(int i)         { lastAtom_ = i;       }
    inline void SetOriginalNum(int i)      { originalResNum_ = i; }
    inline void SetIcode(char c)           { icode_ = c;          }
    inline void SetChainID(char c)         { chainID_ = c;        }
    inline void SetName(NameType const& n) { resname_ = n;        }
    /// \return First atom in residue, indexing from 0
    inline int FirstAtom()        const { return firstAtom_;      }
    /// \return Atom _after_ the last in residue, indexing from 0
    inline int LastAtom()         const { return lastAtom_;       }
    inline int OriginalResNum()   const { return originalResNum_; }
    inline char Icode()           const { return icode_;          }
    inline char ChainID()         const { return chainID_;        }
    inline const char *c_str()    const { return *resname_;       }
    inline NameType const& Name() const { return resname_;        }
    inline int NumAtoms()         const { return (lastAtom_ - firstAtom_); }
    inline bool NameIsSolvent()   const {
      return (resname_=="WAT " || resname_=="HOH " || resname_=="TIP3" ||
              resname_=="SOL ");
    }
    /// Convert 3-letter residue code to single letter.
    static char ConvertResName(std::string const&);
    /// Convert 1-letter residue code to 3 letters.
    static const char* ConvertResName(char);
    /// Convert this residue name to single letter.
    char SingleCharName() const { return ConvertResName( *resname_ ); }
  private:
    NameType resname_;   ///< Residue name.
    int firstAtom_;      ///< Index of first atom (from 0).
    int lastAtom_;       ///< Atom index after last atom in residue.
    int originalResNum_; ///< Original residue number.
    char icode_;         ///< Residue insertion code.
    char chainID_;
};
#endif
