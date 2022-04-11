#ifndef INC_RESIDUE_H
#define INC_RESIDUE_H
#include "NameType.h"
// Class: Residue
/// Hold information for a residue.
class Residue {
  public:
    /// CONSTRUCTOR
    Residue() :
      resname_(""), firstAtom_(0), lastAtom_(0), originalResNum_(0), segID_(-1),
      icode_(' '), chainID_(BLANK_CHAINID_), isTerminal_(false)
    {}
    /// CONSTRUCTOR - Copy given Residue, set first and last atom indices.
    Residue(Residue const& r, int first, int last) :
      resname_(r.resname_), firstAtom_(first), lastAtom_(last),
      originalResNum_(r.originalResNum_), segID_(r.segID_), icode_(r.icode_),
      chainID_(r.chainID_), isTerminal_(false)
    {}
    /// CONSTRUCTOR - Res name, original resnum, icode, chain ID
    Residue(NameType const& n, int r, char ic, char cid) :
      resname_(n), firstAtom_(-1), lastAtom_(-1), originalResNum_(r), segID_(-1),
      icode_(ic), chainID_(cid), isTerminal_(false)
    {}
    /// CONSTRUCTOR - Res name, first atom, last atom, original resnum, icode, chain ID
    Residue(NameType const& n, int first, int last, int r, char ic, char cid) :
      resname_(n), firstAtom_(first), lastAtom_(last),
      originalResNum_(r), segID_(-1), icode_(ic), chainID_(cid), isTerminal_(false)
    {}
    /// CONSTRUCTOR - Res name, original resnum, res icode, segment ID
    Residue(NameType const& n, int r, char i, int s) :
      resname_(n), firstAtom_(-1), lastAtom_(-1), originalResNum_(r), segID_(s),
       icode_(i), chainID_(BLANK_CHAINID_), isTerminal_(false)
    {}
    /// \return True if this residue does not match given residue
    bool operator!=(const Residue& rhs) const {
      return ( originalResNum_ != rhs.originalResNum_ ||
               segID_          != rhs.segID_          ||
               icode_          != rhs.icode_          ||
               chainID_        != rhs.chainID_        ||
               (originalResNum_ == rhs.originalResNum_ && resname_ != rhs.resname_) );
    }
    /// \return Absolute distance in orig. residue numbering
    int AbsResDist(const Residue& rhs) const {
      if (chainID_ != rhs.chainID_) return 99999999; // TODO a real constant
      int dist;
      if (originalResNum_ == rhs.originalResNum_)
      {
        // blank icodes need to be replaced with the one before 'A'
        static const int blankCode = ((int)'A') - 1;
        int code1, code2;
        if (icode_ == ' ')
          code1 = blankCode;
        else
          code1 = (int)icode_;
        if (rhs.icode_ == ' ')
          code2 = blankCode;
        else
          code2 = (int)rhs.icode_;
        dist = code1 - code2;
      } else
        dist = originalResNum_ - rhs.originalResNum_;
      if (dist < 0) dist = -dist;
      return dist;
    }

    inline void SetFirstAtom(int i)        { firstAtom_ = i;      }
    inline void SetLastAtom(int i)         { lastAtom_ = i;       }
    inline void SetOriginalNum(int i)      { originalResNum_ = i; }
    inline void SetSegID(int s)            { segID_ = s;          }
    inline void SetIcode(char c)           { icode_ = c;          }
    inline void SetChainID(char c)         { chainID_ = c;        }
    inline void SetName(NameType const& n) { resname_ = n;        }
    inline void SetTerminal(bool t)        { isTerminal_ = t;     }
    /// \return First atom in residue, indexing from 0
    inline int FirstAtom()        const { return firstAtom_;      }
    /// \return Atom _after_ the last in residue, indexing from 0
    inline int LastAtom()         const { return lastAtom_;       }
    inline int OriginalResNum()   const { return originalResNum_; }
    inline int SegID()            const { return segID_;          }
    inline char Icode()           const { return icode_;          }
    /// \return Chain ID
    inline char ChainId()         const { return chainID_; }
    /// \return True if chain ID is not blank.
    bool HasChainID()             const { return (chainID_ != BLANK_CHAINID_); }
    inline const char *c_str()    const { return *resname_;       }
    inline NameType const& Name() const { return resname_;        }
    inline bool IsTerminal()      const { return isTerminal_;     }
    inline int NumAtoms()         const { return (lastAtom_ - firstAtom_); }
    inline bool NameIsSolvent()   const {
      return (resname_=="WAT " || resname_=="HOH " || resname_=="TIP3" ||
              resname_=="SOL ");
    }
    /// Convert 3-letter residue code to single letter.
    static char ConvertResName(std::string const&);
    /// Convert 1-letter residue code to 3 letters.
    static const char* ConvertResName(char);
    /// \return Default chain ID
    static char DefaultChainID() { return DEFAULT_CHAINID_; }
    /// Convert this residue name to single letter.
    char SingleCharName() const { return ConvertResName( *resname_ ); }
  private:
    /// Character that denotes no chain ID.
    static const char BLANK_CHAINID_;
    /// Chain ID to use if one is desired but no chain ID set.
    static const char DEFAULT_CHAINID_;
    NameType resname_;   ///< Residue name.
    int firstAtom_;      ///< Index of first atom (from 0).
    int lastAtom_;       ///< Atom index after last atom in residue.
    int originalResNum_; ///< Original residue number.
    int segID_;          ///< Segment ID index.
    char icode_;         ///< Residue insertion code.
    char chainID_;       ///< Residue chain ID
    bool isTerminal_;    ///< True if residue was originally a terminal residue
};
#endif
