#ifndef INC_TOPOLOGY_H
#define INC_TOPOLOGY_H
#include <string>
#include <set> // BP_mapType
#include "Atom.h"
#include "Residue.h"
#include "Molecule.h"
#include "ParameterTypes.h"
#include "ParameterSet.h"
#include "Frame.h"
#include "FileName.h"
#include "Range.h"
class AtomMask;
class CharMask;
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
    // ----- Amber Extra Info --------------------
    std::vector<NameType> const& TreeChainClassification() const { return tree_;   }
    std::vector<int>      const& JoinArray()               const { return ijoin_;  }
    std::vector<int>      const& RotateArray()             const { return irotat_; }
    void AllocTreeChainClassification() { tree_.assign(atoms_.size(), "BLA"); }
    void AllocJoinArray()               { ijoin_.assign(atoms_.size(), 0);    }
    void AllocRotateArray()             { irotat_.assign(atoms_.size(), 0);   }
    void SetTreeChainClassification(int idx, NameType const& n) { tree_[idx] = n;   }
    void SetJoinArray(int idx, int j)                           { ijoin_[idx] = j;  }
    void SetRotateArray(int idx, int r)                         { irotat_[idx] = r; }
    // ----- PDB info ----------------------------
    /// Reset all PDB-related info
    void ResetPDBinfo();
    std::vector<char> const& AtomAltLoc()  const { return atom_altloc_;  }
    std::vector<float> const& Occupancy()  const { return occupancy_;    }
    std::vector<float> const& Bfactor()    const { return bfactor_;      }
    std::vector<int> const& PdbSerialNum() const { return pdbSerialNum_; }
    void AllocAtomAltLoc()   { atom_altloc_.assign(atoms_.size(), ' '); }
    void AllocOccupancy()    { occupancy_.assign(atoms_.size(), 1.0);   }
    void AllocBfactor()      { bfactor_.assign(atoms_.size(), 0.0);     }
    void AllocPdbSerialNum() { pdbSerialNum_.assign(atoms_.size(), -1); }
    void SetAtomAltLoc(int idx, char a)  { atom_altloc_[idx] = a;  }
    void SetOccupancy(int idx, float o)  { occupancy_[idx] = o;    }
    void SetBfactor(int idx, float b)    { bfactor_[idx] = b;      }
    void SetPdbSerialNum(int idx, int i) { pdbSerialNum_[idx] = i; }
    void AddAtomAltLoc(char a)  { atom_altloc_.push_back( a );  }
    void AddOccupancy(float o)  { occupancy_.push_back( o );    }
    void AddBfactor(float b)    { bfactor_.push_back( b );      }
    void AddPdbSerialNum(int i) { pdbSerialNum_.push_back( i ); }
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
    /// \return number of residues in the specified molecule.
    int NresInMol(int) const;
    /// Determine molecules based on bond information
    int DetermineMolecules();
    // ----- Bond-specific routines --------------
    size_t Nbonds()                            const { return bonds_.size()+bondsh_.size(); }
    BondArray         const& Bonds()        const { return bonds_;        }
    BondArray         const& BondsH()       const { return bondsh_;       }
    BondParmArray     const& BondParm()     const { return bondparm_;     }
    BondParmType& SetBondParm(int i)              { return bondparm_[i];  }
    void AddBondParm(BondParmType const& b)       { bondparm_.push_back( b ); }
    void AddBond(int i, int j)                    { AddBond(i, j, -1); }
    void AddBond(int, int, int);
    void AddBond(BondType const&, bool);
    void AddBond(int, int, BondParmType const&);
    int RemoveBond(int, int);
    void AssignBondParams(ParmHolder<BondParmType> const&);
    // ----- Angle-specific routines -------------
    size_t Nangles()                           const { return angles_.size()+anglesh_.size(); }
    AngleArray        const& Angles()       const { return angles_;       }
    AngleArray        const& AnglesH()      const { return anglesh_;      }
    AngleParmArray    const& AngleParm()    const { return angleparm_;    }
    AngleParmType& SetAngleParm(int i)            { return angleparm_[i]; }
    void AddAngle(int i, int j, int k)            { AddAngle(i, j, k, -1); }
    void AddAngle(int, int, int, int);
    void AddAngle(AngleType const&, bool);
    void AddAngle(int, int, int, AngleParmType const&); 
    void AssignAngleParams(ParmHolder<AngleParmType> const&);
    // ----- Dihedral-specific routines ----------
    size_t Ndihedrals()                        const { return dihedrals_.size()+dihedralsh_.size(); }
    DihedralArray     const& Dihedrals()    const { return dihedrals_;       }
    DihedralArray     const& DihedralsH()   const { return dihedralsh_;      }
    DihedralParmArray const& DihedralParm() const { return dihedralparm_;    }
    DihedralParmType& SetDihedralParm(int i)      { return dihedralparm_[i]; }
    void AddDihedral(DihedralType const& d)       { AddDihedral(d, -1);      }
    void AddDihedral(DihedralType const&, int);
    void AddDihedral(int i, int j, int k, int l) { AddDihedral(DihedralType(i,j,k,l,-1), -1); }
    void AddDihedral(DihedralType const&, bool);
    void AddDihedral(DihedralType const&, DihedralParmType const&);
    void AssignImproperParams(ParmHolder<DihedralParmType> const&);
    void AssignDihedralParams(DihedralParmHolder const&);
    // ----- CMAP-specific routines --------------
    bool                     HasCmap()      const { return !cmapGrid_.empty(); }
    CmapGridArray     const& CmapGrid()     const { return cmapGrid_;     }
    CmapArray         const& Cmap()         const { return cmap_;         }
    CmapGridType& SetCmapGrid(int idx)            { return cmapGrid_[idx];}
    void AddCmapGrid(CmapGridType const& g) { cmapGrid_.push_back(g); }
    void AddCmapTerm(CmapType const& c)     { cmap_.push_back(c);     }
    // ----- Non-bond routines -------------------
    NonbondParmType  const& Nonbond()        const { return nonbond_;      }
    NonbondParmType&        SetNonbond()           { return nonbond_;      }
    double GetVDWradius(int) const;
    double GetVDWsigma(int) const;
    double GetVDWdepth(int) const;
    /// \return Lennard-Jones 6-12 parameters for given pair of atoms
    inline NonbondType const& GetLJparam(int, int) const;
    void AssignNonbondParams(ParmHolder<AtomType> const&, ParmHolder<NonbondType> const&);
    // ----- Water Cap Info ----------------------
    CapParmType const& Cap()    const { return cap_; }
    CapParmType&       SetCap()       { return cap_; }
    // ----- Amber LES info ----------------------
    LES_ParmType  const& LES()    const { return lesparm_; }
    LES_ParmType&        SetLES()       { return lesparm_; }
    // ----- CHAMBER info ------------------------
    ChamberParmType const& Chamber()        const { return chamber_;      }
    ChamberParmType& SetChamber()                 { return chamber_;      }
    void AddCharmmImproper(DihedralType const&, DihedralParmType const&);
    void AddCharmmImproper(DihedralType const&, int);
    void AddCharmmImproper(DihedralType const& i) { AddCharmmImproper(i, -1); }
    void AssignUBParams(ParmHolder<BondParmType> const&);
    // ----- Misc routines -----------------------
    /// Format: <res name>_<res num>@<atom name>
    std::string TruncResAtomName(int) const;
    /// Format: <res name>@<atom name>
    std::string TruncResNameAtomName(int) const;
    /// Format: <res name>_<res num>@<atom name>_<atom num>
    std::string TruncResAtomNameNum(int) const;
    /// Format: <res name> <res num> <atom name> <atom num>
    std::string ResNameNumAtomNameNum(int) const;
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
    /// \return ParameterSet for this Topology
    ParameterSet GetParameters() const;
    /// Update parameters in this Topology with those in given set.
    int UpdateParams(ParameterSet const&);
    // ----- Print topology info -----------------
    void Summary() const;
    void Brief(const char*) const;
    // ----- Routines to Access/Modify Box info --
    inline Box const& ParmBox()   const { return parmBox_;        }
    void SetParmBox( Box const& bIn )   { parmBox_ = bIn;         }
    void SetBoxFromTraj(Box const&);
    // ----- Setup routines ----------------------
    int AddTopAtom(Atom const&, Residue const&);
    //void StartNewMol();
    /// Perform common final setup: optional molecule determination, renumber residues by molecules
    int CommonSetup(bool, bool);
    /// Perform common final setup with molecule determination on, renumber residues off.
    int CommonSetup() { return CommonSetup(true, false); }
    /// Set up with no residue info TODO deprecate in favor of routine in CommonSetup?
    int Setup_NoResInfo();
    /// Resize for given numbers of atoms/residues etc. Clears any existing data.
    void Resize(Pointers const&);
    // ----- Mask Routines -----------------------
    int SetupIntegerMask(AtomMask &) const;
    int SetupCharMask(CharMask &) const;
    int SetupIntegerMask(AtomMask &, Frame const&) const;
    int SetupCharMask(CharMask &, Frame const&) const;
    /// \return Array of residue numbers selected by given atom mask
    std::vector<int> ResnumsSelectedBy(AtomMask const&) const;
    /// \return Array of molecule numbers selected by given atom mask
    std::vector<int> MolnumsSelectedBy(AtomMask const&) const;
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
    static inline int AddTorsionParm(DihedralParmArray&, DihedralParmType const&);
    bool CheckTorsionRange(DihedralType const& dihIn, const char*) const;
    static inline DihedralType SetTorsionParmIndex(DihedralType const&,
                                                   DihedralParmArray const&,
                                                   int, const char*);
    inline bool CheckExtraSize(size_t, const char*) const;

    void VisitAtom(int, int);
    int RecursiveMolSearch();
    int NonrecursiveMolSearch();
    void ClearMolecules();
    void AtomDistance(int, int, int, std::set<int>&, int) const;
    void DetermineNumExtraPoints();
    int SetSolventInfo();

    int scale_dihedral_K(DihedralArray&, CharMask const&, double, bool);

    Topology* ModifyByMap(std::vector<int> const&, bool) const;
    BondArray StripBondArray(BondArray const&, std::vector<int> const&) const;
    AngleArray StripAngleArray(AngleArray const&, std::vector<int> const&) const;
    DihedralArray StripDihedralArray(DihedralArray const&, std::vector<int> const&) const;
    void StripBondParmArray(BondArray&, std::vector<int>&, BondParmArray&) const;
    void StripBondParmArray(BondArray&, std::vector<int>&, BondParmArray&,
                            BondParmArray const&) const;
    void StripAngleParmArray(AngleArray&, std::vector<int>&, AngleParmArray&) const;
    void StripDihedralParmArray(DihedralArray&, std::vector<int>&, DihedralParmArray&) const;
    void StripDihedralParmArray(DihedralArray&, std::vector<int>&, DihedralParmArray&,
                                DihedralParmArray const&) const;
    inline void AddBondArray(BondArray const&, BondParmArray const&, int);
    inline void AddAngleArray(AngleArray const&, AngleParmArray const&, int);
    inline void AddDihArray(DihedralArray const&, DihedralParmArray const&, int);

    void AssignBondParm(ParmHolder<BondParmType> const&, ParmHolder<int>&, BondArray&, BondParmArray&, const char*);
    void AssignAngleParm(ParmHolder<AngleParmType> const&, ParmHolder<int>&, AngleArray&);
    void AssignImproperParm(ParmHolder<DihedralParmType> const&, ParmHolder<int>&, DihedralArray&);
    void AssignDihedralParm(DihedralParmHolder const&, DihedralArray&);

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
    CmapArray cmap_;                 ///< Hold atom indices and CMAP grid index
    CmapGridArray cmapGrid_;         ///< Hold CMAP grids
    NonbondParmType nonbond_;        ///< Non-bonded parameters
    // Amber-only parameters
    CapParmType cap_;                ///< Water cap information
    LES_ParmType lesparm_;           ///< LES parameters
    ChamberParmType chamber_;        ///< CHAMBER parameters
    // "Extra" Amber atom info
    std::vector<NameType> tree_;     ///< Amber TREE_CHAIN_CLASSIFICATION array
    std::vector<int> ijoin_;         ///< Amber JOIN_ARRAY array
    std::vector<int> irotat_;        ///< Amber IROTAT array
    // PDB info
    std::vector<char> atom_altloc_;  ///< Atom alternate location ID
    std::vector<float> occupancy_;   ///< Atom occupancy
    std::vector<float> bfactor_;     ///< Atom B-factor
    std::vector<int> pdbSerialNum_;  ///< Atom PDB original serial number

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
    Pointers() : natom_(0), nres_(0), nBndParm_(0), nAngParm_(0), nDihParm_(0) {}
    Pointers(int na, int nr, int nbp, int nap, int ndp):
     natom_(na), nres_(nr), nBndParm_(nbp), nAngParm_(nap), nDihParm_(ndp) {} 
  //private:
    int natom_;
    int nres_;
    int nBndParm_;
    int nAngParm_;
    int nDihParm_;
};
#endif
