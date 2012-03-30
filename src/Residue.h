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

    inline const char *c_str() const {
      return *resname_;
    }
    inline void SetName(NameType nameIn) {
      resname_ = nameIn;
    }
    inline NameType Name() const {
      return resname_;
    }
    int OriginalNum();

    void SetFirstAtom(int);
    inline int FirstAtom() const {
      return firstAtom_;
    }

    bool NameIsSolvent();
  private:
    int original_resnum_;
    int firstAtom_;
    NameType resname_;
};
#endif
