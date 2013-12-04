#ifndef INC_TOPOLOGY_H
#define INC_TOPOLOGY_H
#include <string>
#include "Atom.h"
#include "Residue.h"
#include "Molecule.h"
#include "ParameterTypes.h"
#include "AtomMask.h"
#include "Frame.h"
#include "FileName.h"
// Class: Topology
/// Hold information for all atoms
class Topology {
  public:
    Topology();
    // ----- Set internal variables -----
    void SetOffset(double oIn)           { if (oIn > 0.0) offset_ = oIn;  }
    void SetDebug(int dIn)               { debug_ = dIn;                  }
    void SetPindex(int pIn)              { pindex_ = pIn;                 }
    void IncreaseFrames(int fIn)         { nframes_ += fIn;               }
    void SetTag(std::string const& t)    { parmTag_ = t;                  }
    void SetVelInfo(bool v)              { hasVelInfo_ = v;               }
    void SetNrepDim(int n)               { nRepDim_ = n;                  }
    void SetGBradiiSet(std::string const& s) { radius_set_ = s;           }
    void SetParmName(std::string const&, FileName const&);
    void SetReferenceCoords( Frame const& );
    // ----- Return internal variables -----
    std::string const& Tag()       const { return parmTag_;               }
    int Pindex()                   const { return pindex_;                }
    int Natom()                    const { return (int)atoms_.size();     }
    int Nres()                     const { return (int)residues_.size();  }
    int Nmol()                     const { return (int)molecules_.size(); }
    int Nsolvent()                 const { return NsolventMolecules_;     }
    int Nframes()                  const { return nframes_;               }
    int NextraPts()                const { return n_extra_pts_;           }
    bool HasVelInfo()              const { return hasVelInfo_;            }
    int NrepDim()                  const { return nRepDim_;               }
    std::string const& ParmName()         const { return parmName_;       }
    FileName const& OriginalFilename()    const { return fileName_;       }
    std::string const& GBradiiSet()       const { return radius_set_;     }
    bool NoRefCoords()                    const { return (refCoords_.empty()); }
    int FinalSoluteRes() const; // TODO: Replace
    const char *c_str() const;
    // ---- Atom-specific routines -----
    typedef std::vector<Atom>::const_iterator atom_iterator;
    atom_iterator begin()            const { return atoms_.begin(); }
    atom_iterator end()              const { return atoms_.end();   }
    const Atom &operator[](int idx)  const { return atoms_[idx];    }
    std::vector<Atom> const& Atoms() const { return atoms_;         }
    // ----- Residue-specific routines -----
    typedef std::vector<Residue>::const_iterator res_iterator;
    inline res_iterator ResStart() const { return residues_.begin(); }
    inline res_iterator ResEnd()   const { return residues_.end();   }
    const Residue& Res(int idx)    const { return residues_[idx];    }
    // ----- Molecule-specific routines -----
    typedef std::vector<Molecule>::const_iterator mol_iterator;
    inline mol_iterator MolStart() const { return molecules_.begin(); }
    inline mol_iterator MolEnd()   const { return molecules_.end();   }
    const Molecule& Mol(int idx)   const { return molecules_[idx];    }
    // ----- Bond-specific routines -----
    BondArray         const& Bonds()        const { return bonds_;        }
    BondArray         const& BondsH()       const { return bondsh_;       }
    BondParmArray     const& BondParm()     const { return bondparm_;     }
    void AddBond(int,int);
    int SetBondInfo(BondArray const&, BondArray const&, BondParmArray const&);
    // ----- Angle-specific routines -----
    AngleArray        const& Angles()       const { return angles_;       }
    AngleArray        const& AnglesH()      const { return anglesh_;      }
    AngleParmArray    const& AngleParm()    const { return angleparm_;    }
    int SetAngleInfo(AngleArray const&, AngleArray const&, AngleParmArray const&);
    // ----- Dihedral-specific routines -----
    DihedralArray     const& Dihedrals()    const { return dihedrals_;    }
    DihedralArray     const& DihedralsH()   const { return dihedralsh_;   }
    DihedralParmArray const& DihedralParm() const { return dihedralparm_; }
    int SetDihedralInfo(DihedralArray const&, DihedralArray const&, DihedralParmArray const&);
    // ----- Non-bond routines -----
    NonbondParmType   const& Nonbond()      const { return nonbond_;      }
    int SetNonbondInfo(NonbondParmType const&);
    /// \return True if nonbond parameters present.
    bool HasNonbond()                       const { return nonbond_.Ntypes() > 0; }
    /// \return Lennard-Jones 6-12 parameters for given pair of atoms
    inline NonbondType const& GetLJparam(int, int) const;
    // ----- Amber LES info -----
    LES_ParmType      const& LES()          const { return lesparm_;      }
    void SetLES(LES_ParmType const& l)            { lesparm_ = l;         }
    // ----- Amber perturbed parm info -----
    PertParmType      const& Pert()         const { return pert_;         }
    // ----- Amber extra info ----- TODO: Generate automatically, consolidate
    inline const std::vector<double>& Solty()   const { return solty_;  }
    inline const std::vector<NameType>& Itree() const { return itree_;  }
    inline const std::vector<int>& Join()       const { return join_;   }
    inline const std::vector<int>& Irotat()     const { return irotat_; }
    // ----- Misc routines -----
    std::string TruncResAtomName(int) const;
    std::string AtomMaskName(int atom) const;
    std::string TruncResNameNum(int) const;
    int FindAtomInResidue(int, NameType const&) const;
    int FindResidueMaxNatom() const;
    int SoluteAtoms() const;
    int SetSolvent(std::string const&);
    // ----- Print topology info -----
    void Summary() const;
    void Brief(const char*) const;
    void PrintAtomInfo(std::string const&) const;
    void PrintBondInfo(std::string const&) const;
    void PrintAngleInfo(std::string const&) const;
    void PrintDihedralInfo(std::string const&) const;
    void PrintMoleculeInfo(std::string const&) const;
    void PrintResidueInfo(std::string const&) const;
    int PrintChargeMassInfo(std::string const&, int) const;
    // ----- Routines to Access/Modify Box info -----
    inline Box const& ParmBox()   const { return box_;        }
    inline Box::BoxType BoxType() const { return box_.Type(); }
    void SetBox( Box const& bIn )       { box_ = bIn;         }
    // ----- Setup routines -----
    int AddTopAtom(Atom const&, int, NameType const&, const double*);
    void StartNewMol();
    // ----- Amber setup routines -----
    int SetAmberExtra(std::vector<double> const&,std::vector<NameType> const&,
                      std::vector<int> const&,std::vector<int> const&);
    // ----- Common Setup Routines -----
    int CommonSetup(bool);
    // ----- Mask Routines -----
    bool SetupIntegerMask(AtomMask &) const;
    bool SetupCharMask(AtomMask &) const;
    bool SetupIntegerMask(AtomMask &, Frame const&) const;
    bool SetupCharMask(AtomMask &, Frame const&) const;
    // ----- Topology modification routines -----
    void ScaleDihedralK(double);
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
  private:
    void PrintBonds(BondArray const&, AtomMask const&, int&) const;
    void PrintAngles(AngleArray const&, AtomMask const&, int&) const;
    void PrintDihedrals(DihedralArray const&, AtomMask const&, int&) const;
    void SetAtomBondInfo(BondArray const&);
    // NOTE: Use set so that elements are always sorted.
    typedef std::vector< std::set<Atom::AtomicElementType> > BP_mapType;
    void AddBondParam(BondType&, BP_mapType&);
    void AssignBondParameters();
    void GetBondsFromAtomCoords();
    void VisitAtom(int, int);
    int DetermineMolecules();
    void AtomDistance(int, int, int, std::set<int>&) const;
    void DetermineExcludedAtoms();
    void DetermineNumExtraPoints();
    int SetSolventInfo();

    void Mask_SelectDistance( Frame const&, char*, bool, bool, double ) const;
    void Mask_AND(char*,char*) const;
    void Mask_OR(char*,char*) const;
    void Mask_NEG(char*) const;
    void MaskSelectResidues(NameType const&, char *) const;
    void MaskSelectResidues(int, int, char *) const;
    void MaskSelectElements( NameType const&, char* ) const;
    void MaskSelectTypes( NameType const& , char* ) const;
    void MaskSelectAtoms(NameType const&, char*) const;
    void MaskSelectAtoms(int, int, char*) const;
    bool ParseMask(Frame const&, AtomMask &,bool) const;

    Topology* ModifyByMap(std::vector<int> const&, bool) const;
    BondArray StripBondArray(BondArray const&, std::vector<int> const&) const;
    AngleArray StripAngleArray(AngleArray const&, std::vector<int> const&) const;
    DihedralArray StripDihedralArray(DihedralArray const&, std::vector<int> const&) const;

    static const NonbondType LJ_EMPTY;
    std::vector<Atom> atoms_;
    std::vector<Residue> residues_;
    std::vector<Molecule> molecules_;
    FileName fileName_;
    std::string parmTag_;
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
    // Non-bonded parameters
    NonbondParmType nonbond_;
    // Amber-only parameters
    LES_ParmType lesparm_;           ///< LES parameters
    PertParmType pert_;              ///< Atom perturbation info
    // Amber extra info
    std::vector<double> solty_;
    std::vector<NameType> itree_;
    std::vector<int> join_;
    std::vector<int> irotat_;

    Box box_;
    Frame refCoords_;

    double offset_;         ///< Offset used when searching for bonds
    int debug_;
    int NsolventMolecules_;
    int finalSoluteRes_; ///< TODO: Get rid of
    int pindex_;
    int nframes_;
    int n_extra_pts_;
    bool hasVelInfo_; // TODO: This information should be passed separate from Topology
    int nRepDim_;     // TODO: This information should be passed separate from Topology
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
NonbondType const& Topology::GetLJparam(int a1, int a2) const {
  int nbindex = nonbond_.GetLJindex( atoms_[a1].TypeIndex(), atoms_[a2].TypeIndex() );
  if (nbindex < 0) // Means Amber Hbond, return A = B = 0.0
    return LJ_EMPTY;
  return nonbond_.NBarray( nbindex );
}
#endif
