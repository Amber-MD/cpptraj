#ifndef INC_ATOMEXTRA_H
#define INC_ATOMEXTRA_H
#include "NameType.h"
/// Hold PDB anisotropic B factors
class PdbAnisoU {
  public:
    PdbAnisoU() { U_[0] = -1; U_[1] = -1; U_[2] = -1; U_[3] = -1; U_[4] = -1; U_[5] = -1; }
    PdbAnisoU(int u11, int u22, int u33, int u12, int u13, int u23) {
      U_[0] = u11; U_[1] = u22; U_[2] = u33; U_[3] = u12; U_[4] = u13; U_[5] = u23;
    }
    void SetAnisoU(int u11, int u22, int u33, int u12, int u13, int u23) {
      U_[0] = u11; U_[1] = u22; U_[2] = u33; U_[3] = u12; U_[4] = u13; U_[5] = u23;
    }
    int U11() const { return U_[0]; }
    int U22() const { return U_[1]; }
    int U33() const { return U_[2]; }
    int U12() const { return U_[3]; }
    int U13() const { return U_[4]; }
    int U23() const { return U_[5]; }
    bool HasAnisoU() const { return (U_[0] > -1); }
  private:
    int U_[6];
};

/// Hold extra atom info.
class AtomExtra {
  public:
    AtomExtra() : join_(0), irotat_(0), atom_altloc_(' '), occupancy_(1.0), bfactor_(0.0) {}
    AtomExtra(float oc, float bf, char al) :
      join_(0), irotat_(0), atom_altloc_(al), occupancy_(oc), bfactor_(bf) {}
    AtomExtra(NameType const& it, int jo, int ir, char al) :
      itree_(it), join_(jo), irotat_(ir), atom_altloc_(al), occupancy_(1.0), bfactor_(0.0) {}
    // NOTE: This version is for Action_AtomicFluct
    AtomExtra(PdbAnisoU const& a) :
      join_(0), irotat_(0), atom_altloc_(' '), occupancy_(1.0), bfactor_(0.0), anisou_(a) {}

    NameType const& Itree() const { return itree_;       }
    int Join()              const { return join_;        }
    int Irotat()            const { return irotat_;      }
    char AtomAltLoc()       const { return atom_altloc_; }
    float Occupancy()       const { return occupancy_;   }
    float Bfactor()         const { return bfactor_;     }
    PdbAnisoU AnisoU()      const { return anisou_;      }
    void SetAltLoc(char c)           { atom_altloc_ = c;    }
    void SetItree(NameType const& t) { itree_ = t;          }
    void SetJoin(int j)              { join_ = j;           }
    void SetIrotat(int i)            { irotat_ = i;         }
    void SetAnisoU(PdbAnisoU const& a) { anisou_ = a;       }
  private:
    // Amber extra info.
    NameType itree_;
    int join_;
    int irotat_;
    // PDB extra info.
    char atom_altloc_; ///< Alternate location indicator.
    float occupancy_;  ///< Atom occupancy.
    float bfactor_;    ///< Atom B-factor.
    PdbAnisoU anisou_; ///< Anisotropic atom b-factor
};
#endif
