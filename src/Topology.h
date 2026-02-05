#ifndef INC_TOPOLOGY_H
#define INC_TOPOLOGY_H
#include <string>
#include <set> // AtomDistance 
#include "Atom.h"
#include "Residue.h"
#include "Molecule.h"
#include "ParameterTypes.h"
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
    void CopyTopMetadata(Topology const&);
    void SetParmTitle(std::string const&);
    void SetParmName(std::string const&, FileName const&);
    void SetDistMaskRef( Frame const& );
    /// Add forcefield description
    void AddDescription(std::string const& s) { ff_desc_.push_back( s ); }
    // ----- Return internal variables -----------
    int Ipol()                     const { return ipol_;                  }
    int Pindex()                   const { return pindex_;                }
    int Natom()                    const { return (int)atoms_.size();     }
    int Nres()                     const { return (int)residues_.size();  }
    int Nmol()                     const { return (int)molecules_.size(); }
    int Nsolvent()                 const { return NsolventMolecules_;     }
    int NextraPts()                const { return n_extra_pts_;           }
    /// \return Strings describing forcefield
    std::vector<std::string> const& Description() const { return ff_desc_; }

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
    std::vector<Atom>& ModifyAtoms() { return atoms_; }
    Atom& SetAtom(int idx)                             { return atoms_[idx]; }
    /// \return Count of "heavy" atoms (non-hydrogen, non-extra point)
    unsigned int HeavyAtomCount() const;
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
    void AddTreeChainClassification(NameType const& n)          { tree_.push_back(n); }
    void AddJoinArray(int j)                                    { ijoin_.push_back(j); }
    void AddRotateArray(int r)                                  { irotat_.push_back(r); }
    std::vector<NameType>& ModifyTreeChainClassification() { return tree_; }
    std::vector<int>& ModifyJoinArray()                    { return ijoin_; }
    std::vector<int>& ModifyRotateArray()                  { return irotat_; }
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
    std::vector<char>& ModifyAtomAltLoc()  { return atom_altloc_;  }
    std::vector<float>& ModifyOccupancy()  { return occupancy_;    }
    std::vector<float>& ModifyBfactor()    { return bfactor_;      }
    std::vector<int>& ModifyPdbSerialNum() { return pdbSerialNum_; }
    /// Set list of missing residues and residues missing heteroatoms
    void SetMissingResInfo(std::vector<Residue> const&, std::vector<Residue> const&);
    /// \return list of completely missing residues
    std::vector<Residue> const& MissingRes() const { return missingRes_; }
    /// \return list of residues missing heteroatoms
    std::vector<Residue> const& MissingHet() const { return missingHet_; }
    // ----- Residue-specific routines -----------
    typedef std::vector<Residue>::const_iterator res_iterator;
    inline res_iterator ResStart() const { return residues_.begin(); }
    inline res_iterator ResEnd()   const { return residues_.end();   }
    const Residue& Res(int idx)    const { return residues_[idx];    }
    Residue& SetRes(int idx)             { return residues_[idx];    }
    std::vector<Residue> const& Residues() const { return residues_; }
    Range SoluteResidues() const;
    /// Merge residues in given range
    int MergeResidues(int, int);
    // ----- Molecule-specific routines ----------
    typedef std::vector<Molecule>::const_iterator mol_iterator;
    inline mol_iterator MolStart() const { return molecules_.begin(); }
    inline mol_iterator MolEnd()   const { return molecules_.end();   }
    const Molecule& Mol(int idx)   const { return molecules_[idx];    }
    /// \return Array containing solvent molecule numbers
    std::vector<int> SolventMolNums() const;
    /// \return number of residues in the specified molecule.
    int NresInMol(int) const;
    /// Determine molecules based on bond information
    int DetermineMolecules();
    /// Designate all atoms as part of a single molecule.
    int SetSingleMolecule();
    // ----- Bond-specific routines --------------
    size_t Nbonds()                            const { return bonds_.size()+bondsh_.size(); }
    BondArray         const& Bonds()        const { return bonds_;        }
    BondArray         const& BondsH()       const { return bondsh_;       }
    BondParmArray     const& BondParm()     const { return bondparm_;     }
    BondParmType& SetBondParm(int i)              { return bondparm_[i];  }
    BondArray& ModifyBonds() { return bonds_; }
    BondArray& ModifyBondsH() { return bondsh_; }
    BondParmArray& ModifyBondParm() { return bondparm_; }
    void AddBondParm(BondParmType const& b)       { bondparm_.push_back( b ); }
    void AddBond(int i, int j)                    { AddBond(i, j, -1); }
    void AddBond(int, int, int);
    void AddBond(BondType const&, bool);
    void AddBond(int, int, BondParmType const&);
    int RemoveBond(int, int);
    /// Clear bond arrays but not atom connectivity
    void ClearBondArrays();
    // ----- Angle-specific routines -------------
    size_t Nangles()                           const { return angles_.size()+anglesh_.size(); }
    AngleArray        const& Angles()       const { return angles_;       }
    AngleArray        const& AnglesH()      const { return anglesh_;      }
    AngleParmArray    const& AngleParm()    const { return angleparm_;    }
    AngleParmType& SetAngleParm(int i)            { return angleparm_[i]; }
    AngleArray& ModifyAngles() { return angles_; }
    AngleArray& ModifyAnglesH() { return anglesh_; }
    AngleParmArray& ModifyAngleParm() { return angleparm_; }
    void AddAngle(int i, int j, int k)            { AddAngle(i, j, k, -1); }
    void AddAngle(int, int, int, int);
    void AddAngle(AngleType const&, bool);
    void AddAngle(int, int, int, AngleParmType const&);
    // ----- Dihedral-specific routines ----------
    size_t Ndihedrals()                        const { return dihedrals_.size()+dihedralsh_.size(); }
    DihedralArray     const& Dihedrals()    const { return dihedrals_;       }
    DihedralArray     const& DihedralsH()   const { return dihedralsh_;      }
    DihedralParmArray const& DihedralParm() const { return dihedralparm_;    }
    DihedralParmType& SetDihedralParm(int i)      { return dihedralparm_[i]; }
    DihedralArray& ModifyDihedrals() { return dihedrals_; }
    DihedralArray& ModifyDihedralsH() { return dihedralsh_; }
    DihedralParmArray& ModifyDihedralParm() { return dihedralparm_; }
    void AddDihedral(DihedralType const& d)       { AddDihedral(d, -1);      }
    void AddDihedral(DihedralType const&, int);
    void AddDihedral(int i, int j, int k, int l) { AddDihedral(DihedralType(i,j,k,l,-1), -1); }
    void AddDihedral(DihedralType const&, bool);
    void AddDihedral(DihedralType const&, DihedralParmType const&);
    /// \return Index of existing/added dihedral parameter in given array
    static int addTorsionParm(DihedralParmArray&, DihedralParmType const&);

    // ----- CMAP-specific routines --------------
    bool                     HasCmap()      const { return !cmapGrid_.empty(); }
    CmapGridArray     const& CmapGrid()     const { return cmapGrid_;     }
    CmapArray         const& Cmap()         const { return cmap_;         }
    CmapGridType& SetCmapGrid(int idx)            { return cmapGrid_[idx];}
    CmapGridArray& ModifyCmapGrid() { return cmapGrid_; }
    CmapArray& ModifyCmap() { return cmap_; }
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
    /// \return Lennard-Jones 6-12 parameters for given pair of atoms, set the LJC parameter
    inline NonbondType const& GetLJCparam(double&, int, int) const;
    /// \return True if any charge is non-zero
    bool HasChargeInfo() const;
    /// \return Total charge
    double TotalCharge() const;
    /// Redistribute charge on atoms in topology to match target total charge
    int RedistributeCharge(double);
    // ----- Water Cap Info ----------------------
    CapParmType const& Cap()    const { return cap_; }
    CapParmType&       SetCap()       { return cap_; }
    // ----- Amber LES info ----------------------
    LES_ParmType  const& LES()    const { return lesparm_; }
    LES_ParmType&        SetLES()       { return lesparm_; }
    // ----- CHAMBER info ------------------------
    /// \return true if any CHARMM parameters are set (based on indices).
    bool HasChamber() const;

    /// Reserve space for future UB terms
    void ReserveUBterms(unsigned int n)       { ub_.reserve( n );                 }
    /// Resize for given number of UB parameters
    void ResizeUBparm(unsigned int n)         { ubparm_.resize( n );              }
    /// Add a Urey-Bradley term
    void AddUBterm(BondType const& bnd)       { ub_.push_back( bnd );             }
    /// Set specified Urey-Bradley parameter
    BondParmType& SetUBparm(unsigned int idx) { return ubparm_[idx];              }
    /// \return Array of Urey-Bradley terms
    BondArray         const& UB()       const { return ub_;           }
    /// \return Array of Urey-Bradley parameters
    BondParmArray     const& UBparm()   const { return ubparm_;       }
    /// \return Modifiable array of UB terms
    BondArray& ModifyUB() { return ub_; }
    /// \return Modifiable array of UB parameters
    BondParmArray& ModifyUBparm() { return ubparm_; }

    /// Reserve space for future improper terms
    void ReserveImproperTerms(unsigned int n)         { impropers_.reserve( n );     }
    /// Resize for given number of improper parameters
    void ResizeImproperParm(unsigned int n)           { improperparm_.resize( n );   }
    /// \return Modifiable improper array
    DihedralArray& ModifyImpropers() { return impropers_; }
    /// \return Modifiable improper parm array
    DihedralParmArray& ModifyImproperParm() { return improperparm_; }
    /// Add improper term
    void AddImproperTerm(DihedralType const& dih)     { impropers_.push_back( dih ); }
    /// Set specified improper parameter
    DihedralParmType& SetImproperParm(unsigned int i) { return improperparm_[i];     }
    /// \return Array of impropers
    DihedralArray     const& Impropers()        const { return impropers_;    }
    /// \return Array of improper parameters
    DihedralParmArray const& ImproperParm()     const { return improperparm_; }


    void AddCharmmImproper(DihedralType const&, DihedralParmType const&);
    void AddCharmmImproper(DihedralType const&, int);
    void AddCharmmImproper(DihedralType const& i) { AddCharmmImproper(i, -1); }
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
    /// LEaP Format: R<res #>.A<at #>
    std::string LeapName(int) const;
    /// Format: <res name>:<res num> 
    std::string TruncResNameNum(int) const;
    /// Format: <res name>_<onum>[_<id>]
    std::string TruncResNameOnumId(int) const;
    /// Format: <res name>_<onum>[_<id>]@<atom name>
    std::string TruncAtomResNameOnumId(int) const;
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
    /// Add an atom to the topology inside given residue
    int AddTopAtom(Atom const&, Residue const&);
    /// Add specified residues from given topology; optionally make them solvent
    int AddResidues(Topology const&, std::vector<int> const&, Frame&, Frame const&, bool);
    /// Add an atom to the topology inside given residue and molecule number
    int addTopAtom(Atom const&, Residue const&, unsigned int, bool);
    //void StartNewMol();
    /// Perform common final setup: optional molecule determination, renumber residues by molecules
    int CommonSetup(bool, bool);
    /// Perform common final setup with molecule determination on, renumber residues off
    int CommonSetup() { return CommonSetup(true, false); }
    /// Set up with no residue info TODO deprecate in favor of routine in CommonSetup?
    int Setup_NoResInfo();
    /// Resize for given numbers of atoms/residues etc. Clears any existing data.
    void Resize(Pointers const&);
    /// Count the number of extra points in the Topology TODO just have it a part of AddTopAtom
    void DetermineNumExtraPoints();
    // ----- Mask Routines -----------------------
    int SetupIntegerMask(AtomMask &) const;
    int SetupCharMask(CharMask &) const;
    int SetupIntegerMask(AtomMask &, Frame const&) const;
    int SetupCharMask(CharMask &, Frame const&) const;
    /// \return Array of residue numbers selected by given atom mask
    std::vector<int> ResnumsSelectedBy(AtomMask const&) const;
    /// \return Array of molecule numbers selected by given atom mask
    std::vector<int> MolnumsSelectedBy(AtomMask const&) const;
    /// \return True if the total mass of selected atoms is zero
    bool MaskHasZeroMass(AtomMask const&) const;
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
    /// Append topology to this one. Takes verbosity, reduce bond params, reduce angle params
    int AppendTop(Topology const&, int, bool, bool);
    /// Append topology to this one.  No verbosity, no reduce bond/angle params. PYTRAJ relies on this.
    int AppendTop(Topology const&);
    /// Split selected atoms in a residue into a new residue, populate the atom map
    int SplitResidue(AtomMask const&, NameType const&, std::vector<int>&);
    /// Split selected atoms in a residue into a new residue
    int SplitResidue(AtomMask const&, NameType const&);
    
  private:
    typedef std::vector<std::string> Sarray;

    /// \return Index of existing/added bond parameter in given array
    static inline int addBondParm(BondParmArray&, BondParmType const&);
    /// \return Index of existing/added angle parameter in given array
    static inline int addAngleParm(AngleParmArray&, AngleParmType const&);
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
    // Charmm parameters
    Sarray ff_desc_;                 ///< Force field descriptions
    BondArray ub_;                   ///< Urey-Bradley terms
    BondParmArray ubparm_;           ///< Urey-Bradley parameters
    DihedralArray impropers_;        ///< Improper terms
    DihedralParmArray improperparm_; ///< Improper parameters
    // "Extra" Amber atom info
    std::vector<NameType> tree_;     ///< Amber TREE_CHAIN_CLASSIFICATION array
    std::vector<int> ijoin_;         ///< Amber JOIN_ARRAY array
    std::vector<int> irotat_;        ///< Amber IROTAT array
    // PDB info
    std::vector<char> atom_altloc_;  ///< Atom alternate location ID
    std::vector<float> occupancy_;   ///< Atom occupancy
    std::vector<float> bfactor_;     ///< Atom B-factor
    std::vector<int> pdbSerialNum_;  ///< Atom PDB original serial number
    std::vector<Residue> missingRes_; ///< List of residues missing from PDB
    std::vector<Residue> missingHet_; ///< List of residues missing heteroatoms in PDB

    Box parmBox_;
    Frame refCoords_;       ///< Internal reference coords for distance-based masks

    int debug_;
    int ipol_;              ///< 0 if fixed charge, 1 if polarizable
    int NsolventMolecules_; ///< Number of molecules marked SOLVENT
    int pindex_;            ///< Internal index used to ID Topology 
    int n_extra_pts_;       ///< Number of extra points.
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
NonbondType const& Topology::GetLJparam(int a1, int a2) const {
  int nbindex = nonbond_.GetLJindex( atoms_[a1].TypeIndex(), atoms_[a2].TypeIndex() );
  if (nbindex < 0) // Means Amber Hbond, return A = B = 0.0
    return LJ_EMPTY;
  return nonbond_.NBarray( nbindex );
}
NonbondType const& Topology::GetLJCparam(double& LJC, int a1, int a2) const {
  int nbindex = nonbond_.GetLJindex( atoms_[a1].TypeIndex(), atoms_[a2].TypeIndex() );
  if (nbindex < 0) // Means Amber Hbond, return A = B = 0.0
    return LJ_EMPTY;
  if (nonbond_.Has_C_Coeff())
    LJC = nonbond_.LJC_Array(nbindex);
  else
    LJC = 0;
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
