#ifndef INC_RESIDUE_H
#define INC_RESIDUE_H
#include "NameType.h"
// Class: Residue
/// Hold information for a residue.
class Residue {
  public:
    Residue();
    Residue(int, NameType);
    Residue(NameType, int);
    Residue(int);

    bool operator==(const Residue &);
    bool operator!=(const Residue &);

    inline int OriginalNum()         { return original_resnum_; }
    inline int FirstAtom() const     { return firstAtom_;       }
    inline const char *c_str() const { return *resname_;        }
    inline NameType Name() const     { return resname_;         }

    void SetFirstAtom(int);
    inline void SetName(NameType nameIn) { resname_ = nameIn; }

    bool NameIsSolvent();
  private:
    int original_resnum_;
    int firstAtom_;
    NameType resname_;
};
#endif
