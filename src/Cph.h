#ifndef INC_CPH_H
#define INC_CPH_H
#include <vector>
#include "NameType.h"
namespace Cph {

/// Indicate full record written
const int FULL_RECORD = -2;

/// Indicate partial record written
const int PARTIAL_RECORD = -1;

/// Hold information for a titratable residue.
class CpRes {
    typedef std::vector<int> Iarray;
    typedef std::vector<bool> Barray;
  public:
    CpRes() : resid_(-1) {}
    /// Res name, num, protcnts, max prot
    CpRes(NameType const&, int, Iarray const&, int);
    /// COPY
    CpRes(CpRes const&);
    /// ASSIGN
    CpRes& operator=(CpRes const&);
    //void Print() const;
    bool IsProtonated(int s)  const { return protonated_[s]; }
    int Nprotons(int s)       const { return protcnts_[s];   }
    NameType const& Name()    const { return resname_;       }
    int Num()                 const { return resid_;         }
    /// Indicates full record written.
    static const int FULL_RECORD;
    /// Indicate partial record written.
    static const int PARTIAL_RECORD;
  private:
    NameType resname_;  ///< Residue name.
    int resid_;         ///< Residue number.
    Iarray protcnts_;   ///< Hold protonation count for each state.
    Barray protonated_; ///< True if state is protonated.
};
} // END namepsace Cph
#endif
