#ifndef INC_TOPOLOGY_H
#define INC_TOPOLOGY_H
#include <string>
#include "Atom.h"
#include "AtomExtra.h"
#include "Residue.h"
#include "Molecule.h"
#include "ParameterTypes.h"
#include "AtomMask.h"
#include "CharMask.h"
#include "Frame.h"
#include "FileName.h"
#include "Range.h"
/// Hold information for all atoms
class Topology {
  public:
    Topology();
    /// Can be used with Resize() to reserve Topology memory.
    class Pointers;
    // ----- Set internal variables --------------
    void SetDebug(int dIn)                   { debug_ = dIn;                 }
    void SetIpol(int iIn)                    { ipol_ = iIn;                  }
    void SetPindex(int pIn)                  { pindex_ = pIn;                }
    void SetGBradiiSet(std::string const& s) { radius_set_ = s;              }
    void SetParmName(std::string const&, FileName const&);
    void SetDistMaskRef( Frame const& );
    /// Set value of NATYP from Amber Topology. Only needed for Amber.
    void SetNatyp(int n)                     { n_atom_types_ = n;            }
    // ----- Return internal variables -----------
    int Ipol()                     const { return ipol_;                  }
    int Pindex()                   const { return pindex_;                }
    int Natom()                    const { return (int)atoms_.size();     }
    int Nres()                     const { return (int)residues_.size();  }
    int Nmol()                     const { return (int)molecules_.size(); }
    int Nsolvent()                 const { return NsolventMolecules_;     }
    int NextraPts()                const { return n_extra_pts_;           }
    inline int NatomTypes()        const { return n_atom_types_;          }
    std::string const& ParmName()         const { return parmName_;       }
    FileName const& OriginalFilename()    const { return fileName_;       }
    std::string const& GBradiiSet()       const { return radius_set_;     }
    const char *c_str() const; //FIXME rename
    // ---- Atom-specific routines ---------------
    typedef std::vector<Atom>::const_iterator atom_iterator;
    atom_iterator begin()                        const { return atoms_.begin(); }
    atom_iterator end()                          const { return atoms_.end();   }
    const Atom &operator[](int idx)              const { return atoms_[idx];    }
    std::vector<Atom> const& Atoms()             const { return atoms_;         }
    Atom& SetAtom(int idx)                             { return atoms_[idx]; }
    typedef std::vector<AtomExtra>::const_iterator extra_iterator;
    extra_iterator extraBegin()                  const { return extra_.begin(); }
    extra_iterator extraEnd()                    const { return extra_.end();   }
    inline const std::vector<AtomExtra>& Extra() const { return extra_;         }
    AtomExtra& SetExtraAtomInfo(int idx)               { return extra_[idx];    }
    // ----- Residue-specific routines -----------
    typedef std::vector<Residue>::const_iterator res_iterator;
    inline res_iterator ResStart() const { return residues_.begin(); }
    inline res_iterator ResEnd()   const { return residues_.end();   }
    const Residue& Res(int idx)    const { return residues_[idx];    }
    Residue& SetRes(int idx)             { return residues_[idx];    }
    Range SoluteResidues() const;
    // ----- Molecule-specific routines ----------
    typedef std::vector<Molecule>::const_iterator mol_iterator;
    inline mol_iterator MolStart() const { return molecules_.begin(); }
    inline mol_iterator MolEnd()   const { return molecules_.end();   }
    const Molecule& Mol(int idx)   const { return molecules_[idx];    }
    // ----- Bond-specific routines --------------
    BondArray         const& Bonds()        const { return bonds_;        }
    BondArray         const& BondsH()       const { return bondsh_;       }
    BondParmArray     const& BondParm()     const { return bondparm_;     }
    BondParmType& SetBondParm(int i)              { return bondparm_[i];  }
    void AddBondParm(BondParmType const& b)       { bondparm_.push_back( b ); }
    void AddBond(int i, int j)                    { AddBond(i, j, -1); }
    void AddBond(int, int, int);
    void AddBond(BondType const&, bool);
    // ----- Angle-specific routines -------------
    AngleArray        const& Angles()       const { return angles_;       }
    AngleArray        const& AnglesH()      const { return anglesh_;      }
    AngleParmArray    const& AngleParm()    const { return angleparm_;    }
    AngleParmType& SetAngleParm(int i)            { return angleparm_[i]; }
    void AddAngle(int, int, int);
    void AddAngle(AngleType const&, bool);
    // ----- Dihedral-specific routines ----------
    DihedralArray     const& Dihedrals()    const { return dihedrals_;       }
    DihedralArray     const& DihedralsH()   const { return dihedralsh_;      }
    DihedralParmArray const& DihedralParm() const { return dihedralparm_;    }
    DihedralParmType& SetDihedralParm(int i)      { return dihedralparm_[i]; }
    void AddDihedral(int, int, int, int);
    void AddDihedral(DihedralType const&, bool);
    // ----- Non-bond routines -------------------
    NonbondParmType  const& Nonbond()        const { return nonbond_;      }
    NonbondParmType&        SetNonbond()           { return nonbond_;      }
    double GetVDWradius(int) const;
    double GetVDWdepth(int) const;
    /// \return Lennard-Jones 6-12 parameters for given pair of atoms
    inline NonbondType const& GetLJparam(int, int) const;
    // ----- Water Cap Info ----------------------
    CapParmType const& Cap()    const { return cap_; }
    CapParmType&       SetCap()       { return cap_; }
    // ----- Amber LES info ----------------------
    LES_ParmType  const& LES()    const { return lesparm_; }
    LES_ParmType&        SetLES()       { return lesparm_; }
    // ----- CHAMBER info ------------------------
    ChamberParmType const& Chamber()        const { return chamber_;      }
    ChamberParmType& SetChamber()                 { return chamber_;      }
    // ----- Misc routines -----------------------
    /// Format: <res name><res num>@<atom name>
    std::string TruncResAtomName(int) const;
    /// Format: :<res num>@<atom name>
    std::string AtomMaskName(int) const;
    /// Format: <atom name>_<atom num>
    std::string TruncAtomNameNum(int) const;
    /// Format: <res name>:<res num> 
    std::string TruncResNameNum(int) const;
    /// \return index of atom with given name in specified residue.
    int FindAtomInResidue(int, NameType const&) const;
    /// Mark all molecules matching given mask expression as solvent.
    int SetSolvent(std::string const&);
    // ----- Print topology info -----------------
    void Summary() const;
    void Brief(const char*) const;
    // ----- Routines to Access/Modify Box info --
    inline Box const& ParmBox()   const { return parmBox_;        }
    void SetParmBox( Box const& bIn )   { parmBox_ = bIn;         }
    void SetBoxFromTraj(Box const&);
    // ----- Setup routines ----------------------
    int AddTopAtom(Atom const&, Residue const&);
    void AddExtraAtomInfo(AtomExtra const& ex) { extra_.push_back(ex); } // FIXME bounds check
    void StartNewMol();
    int CommonSetup(bool);
    int CommonSetup() { return CommonSetup(true); }
    void ResetPDBinfo();
    int Setup_NoResInfo();
    /// Resize for given numbers of atoms/residues etc. Clears any existing data.
    void Resize(Pointers const&);
    // ----- Mask Routines -----------------------
    int SetupIntegerMask(AtomMask &) const;
    int SetupCharMask(CharMask &) const;
    int SetupIntegerMask(AtomMask &, Frame const&) const;
    int SetupCharMask(CharMask &, Frame const&) const;
    // ----- Topology modification routines ------
    int ScaleDihedralK(double, std::string const&, bool);
    /// Strip atoms outside given mask, do not keep parameters.
    Topology* partialModifyStateByMask(AtomMask const& m) const {
      return ModifyByMap(m.Selected(), false);
    }
    /// Strip atoms outside given mask.
    Topology* modifyStateByMask(AtomMask const& m) const {
      return ModifyByMap(m.Selected(), true);
    }
    /// Rearrange atoms from given map, Map[newatom]=oldatom
    Topology* ModifyByMap(std::vector<int> const& m) const {
      return ModifyByMap(m, true);
    }
    /// Append topology to this one.
    int AppendTop( Topology const& );
  private:
    void SetAtomBondInfo(BondArray const&);
    // NOTE: Use set so that elements are always sorted.
    typedef std::vector< std::set<Atom::AtomicElementType> > BP_mapType;
    void AddBondParam(BondType&, BP_mapType&);
    void AssignBondParameters();
    void VisitAtom(int, int);
    int RecursiveMolSearch();
    int NonrecursiveMolSearch();
    void ClearMolecules();
    int DetermineMolecules();
    void AtomDistance(int, int, int, std::set<int>&) const;
    void DetermineExcludedAtoms();
    void DetermineNumExtraPoints();
    int SetSolventInfo();

    int scale_dihedral_K(DihedralArray&, CharMask const&, double, bool);

    Topology* ModifyByMap(std::vector<int> const&, bool) const;
    BondArray StripBondArray(BondArray const&, std::vector<int> const&) const;
    AngleArray StripAngleArray(AngleArray const&, std::vector<int> const&) const;
    DihedralArray StripDihedralArray(DihedralArray const&, std::vector<int> const&) const;
    void StripBondParmArray(BondArray&, std::vector<int>&, BondParmArray&) const;
    void StripAngleParmArray(AngleArray&, std::vector<int>&, AngleParmArray&) const;
    void StripDihedralParmArray(DihedralArray&, std::vector<int>&, DihedralParmArray&) const;
    inline void AddBondArray(BondArray const&, int);

    static const NonbondType LJ_EMPTY;
    std::vector<Atom> atoms_;
    std::vector<Residue> residues_;
    std::vector<Molecule> molecules_;
    // NOTE: Filename is stored to enable things like 'strip outprefix'
    FileName fileName_; 
    std::string parmName_;
    std::string radius_set_;

    BondArray bonds_;
    BondArray bondsh_;
    BondParmArray bondparm_;
    AngleArray angles_;
    AngleArray anglesh_;
    AngleParmArray angleparm_;
    DihedralArray dihedrals_;
    DihedralArray dihedralsh_;
    DihedralParmArray dihedralparm_;
    NonbondParmType nonbond_;        ///< Non-bonded parameters
    // Amber-only parameters
    CapParmType cap_;                ///< Water cap information
    LES_ParmType lesparm_;           ///< LES parameters
    ChamberParmType chamber_;        ///< CHAMBER parameters
    // Extra atom info
    std::vector<AtomExtra> extra_;

    Box parmBox_;
    Frame refCoords_;       ///< Internal reference coords for distance-based masks

    int debug_;
    int ipol_;              ///< 0 if fixed charge, 1 if polarizable
    int NsolventMolecules_; ///< Number of molecules marked SOLVENT
    int pindex_;            ///< Internal index used to ID Topology 
    int n_extra_pts_;       ///< Number of extra points.
    int n_atom_types_;      ///< Number of unique atom types.
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
NonbondType const& Topology::GetLJparam(int a1, int a2) const {
  int nbindex = nonbond_.GetLJindex( atoms_[a1].TypeIndex(), atoms_[a2].TypeIndex() );
  if (nbindex < 0) // Means Amber Hbond, return A = B = 0.0
    return LJ_EMPTY;
  return nonbond_.NBarray( nbindex );
}
// -----------------------------------------------------------------------------
class Topology::Pointers {
  public:
    Pointers() : natom_(0), nres_(0), nextra_(0), nBndParm_(0), nAngParm_(0),
                 nDihParm_(0) {}
    Pointers(int na, int nr, int ne, int nbp, int nap, int ndp):
     natom_(na), nres_(nr), nextra_(ne), nBndParm_(nbp), nAngParm_(nap),
     nDihParm_(ndp) {} 
  //private:
    int natom_;
    int nres_;
    int nextra_;
    int nBndParm_;
    int nAngParm_;
    int nDihParm_;
};
#endif
