#ifndef INC_CPHRESIDUE_H
#define INC_CPHRESIDUE_H
#include <vector>
#include "NameType.h"
/// Hold information for a titratable residue.
class CphResidue {
    typedef std::vector<int> Iarray;
    typedef std::vector<bool> Barray;
  public:
    CphResidue() : resid_(-1) {}
    /// Res name, num, protcnts, max prot
    CphResidue(NameType const&, int, Iarray const&, int);
    /// COPY
    CphResidue(CphResidue const&);
    /// ASSIGN
    CphResidue& operator=(CphResidue const&);
    void Print() const;
    bool IsProtonated(int s)  const { return protonated_[s]; }
    int Nprotons(int s)       const { return protcnts_[s];   }
    NameType const& Name()    const { return resname_;       }
    int Num()                 const { return resid_;         }
  private:
    NameType resname_;  ///< Residue name.
    int resid_;         ///< Residue number.
    Iarray protcnts_;   ///< Hold protonation count for each state.
    Barray protonated_; ///< True if state is protonated.
};
#endif
