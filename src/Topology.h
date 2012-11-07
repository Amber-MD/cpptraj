#ifndef INC_TOPOLOGY_H
#define INC_TOPOLOGY_H
#include <string>
#include "Atom.h"
#include "Residue.h"
#include "Molecule.h"
#include "AtomMask.h"
#include "Box.h"
#include "Frame.h"
// Class: Topology
/// Hold information for all atoms
class Topology {
  public:
    Topology();
    // ----- Set internal variables -----
    void SetOffset(double oIn);
    void SetDebug(int);
    void SetParmName(std::string const&, std::string const&);
    void SetGBradiiSet(std::string&);
    void SetPindex(int);
    void SetReferenceCoords( Frame* ); // TODO: Pass in frame reference
    // ----- Return internal variables -----
    int FinalSoluteRes() const;
    int Pindex() const;
    int Natom() const;
    int Nres() const;
    int Nmol() const;
    int Nsolvent() const;
    int Nframes() const;
    int Ntypes() const;
    void IncreaseFrames(int);
    const char *c_str();
    std::string ParmName() const;
    std::string OriginalFilename();
    std::string GBradiiSet() const;
    // ---- Atom-specific routines -----
    typedef std::vector<Atom>::iterator iterator;
    inline iterator begin() { return atoms_.begin(); }
    inline iterator end()   { return atoms_.end();   }
    typedef std::vector<Atom>::const_iterator atom_iterator;
    inline atom_iterator begin() const { return atoms_.begin(); }
    inline atom_iterator end()   const { return atoms_.end();   }
    atom_iterator ResAtomStart(int) const;
    atom_iterator ResAtomEnd(int) const;
    atom_iterator MolAtomStart(int) const;
    atom_iterator MolAtomEnd(int) const;
    const Atom &operator[](int) const;
    std::vector<Atom> const& Atoms() { return atoms_; }
    // ----- Residue-specific routines -----
    typedef std::vector<Residue>::const_iterator res_iterator;
    inline res_iterator ResStart() const { return residues_.begin(); }
    inline res_iterator ResEnd() const   { return residues_.end();   }
    const Residue& Res(int idx)          { return residues_[idx];    }
    int ResFirstAtom(int);
    int ResLastAtom(int);
    int ResSize(int); 
    // ----- Molecule-specific routines -----
    typedef std::vector<Molecule>::const_iterator mol_iterator;
    inline mol_iterator MolStart() const { return molecules_.begin(); }
    inline mol_iterator MolEnd() const   { return molecules_.end();   }
    const Molecule& Mol(int idx)         { return molecules_[idx];    }
    // ----- Bond-specific routines -----
    inline const std::vector<int>& Bonds() const { return bonds_; }
    inline const std::vector<int>& BondsH() const { return bondsh_; }
    inline const std::vector<double>& BondRk() const { return bondrk_; }
    inline const std::vector<double>& BondReq() const { return bondreq_; }
    int GetBondParamIdx( int, double &, double &);
    double GetBondedCutoff(int, int);
    // ----- Angle-specific routines -----
    inline const std::vector<int>& Angles() const { return angles_; }
    inline const std::vector<int>& AnglesH() const { return anglesh_; }
    inline const std::vector<double>& AngleTk() const { return angletk_; }
    inline const std::vector<double>& AngleTeq() const { return angleteq_; }
    // ----- Dihedral-specific routines -----
    inline const std::vector<int>& Dihedrals() const { return dihedrals_; }
    inline const std::vector<int>& DihedralsH() const { return dihedralsh_; }
    inline const std::vector<double>& DihedralPk() const { return dihedralpk_; }
    inline const std::vector<double>& DihedralPn() const { return dihedralpn_; }
    inline const std::vector<double>& DihedralPhase() const { return dihedralphase_; }
    inline const std::vector<double>& SCEE() const { return scee_; }
    inline const std::vector<double>& SCNB() const { return scnb_; }
    // ----- Amber Hbond info -----
    inline const std::vector<double>& Asol() const { return asol_; }
    inline const std::vector<double>& Bsol() const { return bsol_; }
    inline const std::vector<double>& HBcut() const { return hbcut_; }
    // ----- Amber extra info ----- TODO: Generate automatically
    inline const std::vector<double>& Solty() const { return solty_; }
    inline const std::vector<NameType>& Itree() const { return itree_; }
    inline const std::vector<int>& Join() const { return join_; }
    inline const std::vector<int>& Irotat() const { return irotat_; }
    // ----- Non-bond routines -----
    inline const std::vector<int>& NB_index() const { return nbindex_; }
    inline const std::vector<double>& LJA() const { return lja_; }
    inline const std::vector<double>& LJB() const { return ljb_; }
    // ----- Misc routines -----
    int ResAtomRange(int, int *, int *);
    const char* ResidueName(int); // TODO: Make obsolete
    std::string ResAtomName(int);
    std::string TruncResAtomName(int);
    std::string ResNameNum(int);
    std::string TruncResNameNum(int);
    int FindAtomInResidue(int, NameType);
    int FindResidueMaxNatom() const;
    int SoluteAtoms();
    //double *Mass();
    int SetSolvent(std::string const&);
    // ----- Print topology info -----
    void Summary();
    void ParmInfo();
    void PrintAtomInfo(std::string const&);
    void PrintBondInfo();
    void PrintMoleculeInfo();
    void PrintResidueInfo();
    // ----- Routines to Access/Modify Box info -----
    inline Box const& ParmBox() const { return box_;        }
    inline Box::BoxType BoxType()     { return box_.Type(); }
    void SetBox( Box const& bIn )     { box_ = bIn;         }
    // ----- PDB/Mol2 etc setup routines -----
    void AddAtom(Atom, Residue, const double*);
    void StartNewMol();
    // ----- Amber setup routines -----
    int CreateAtomArray(std::vector<NameType>&, std::vector<double>&,
                        std::vector<int>&, std::vector<double>&,
                        std::vector<int>&, std::vector<NameType>&,
                        std::vector<double>&,std::vector<double>&, 
                        std::vector<NameType>&, std::vector<int>&);
    int SetBondInfo(std::vector<int> &, std::vector<int> &,
                    std::vector<double>&,std::vector<double>&);
    int SetAngleInfo(std::vector<int> &, std::vector<int> &,
                    std::vector<double>&,std::vector<double>&);
    int SetDihedralInfo(std::vector<int> &, std::vector<int> &,
                    std::vector<double>&,std::vector<double>&,
                    std::vector<double>&,
                    std::vector<double>&,std::vector<double>&);
    int SetAmberHbond(std::vector<double>&,std::vector<double>&,std::vector<double>&);
    int SetAmberExtra(std::vector<double>&,std::vector<NameType> &,
                      std::vector<int> &,std::vector<int> &);
    int SetNonbondInfo(int, std::vector<int>& nbindex, 
                       std::vector<double>&, std::vector<double>&);
    // ----- Common Setup Routines -----
    int CommonSetup(bool);
    void ClearBondInfo();
    void AddBond(int,int);
    // ----- Mask Routines -----
    bool SetupIntegerMask(AtomMask &);
    bool SetupCharMask(AtomMask &);
    bool SetupIntegerMask(AtomMask &, Frame &);
    bool SetupCharMask(AtomMask &, Frame &);
    // ----- Topology modification routines -----
    Topology *modifyStateByMask(AtomMask &);
    Topology *ModifyByMap(std::vector<int>&);

  private:
    std::vector<Atom> atoms_;
    std::vector<Residue> residues_;
    std::vector<Molecule> molecules_;
    std::string fileName_;
    std::string parmName_;
    std::string radius_set_;

    std::vector<int> bonds_;
    std::vector<int> bondsh_;
    std::vector<double> bondrk_;
    std::vector<double> bondreq_;

    std::vector<int> angles_;
    std::vector<int> anglesh_;
    std::vector<double> angletk_;
    std::vector<double> angleteq_;

    std::vector<int> dihedrals_;
    std::vector<int> dihedralsh_;
    std::vector<double> dihedralpk_;
    std::vector<double> dihedralpn_;
    std::vector<double> dihedralphase_;
    std::vector<double> scee_;
    std::vector<double> scnb_;

    std::vector<double> asol_;
    std::vector<double> bsol_;
    std::vector<double> hbcut_;

    std::vector<double> solty_;
    std::vector<NameType> itree_;
    std::vector<int> join_;
    std::vector<int> irotat_;

    std::vector<int> nbindex_;
    std::vector<double> lja_;
    std::vector<double> ljb_;

    Box box_;
    Frame refCoords_;

    double offset_;         ///< Offset used when searching for bonds
    int debug_;
    int topology_error_;
    int NsolventMolecules_;
    int finalSoluteRes_;
    int pindex_;
    int nframes_;
    int ntypes_; // This is stored for the purpose of checking array sizes

    void PrintBonds(std::vector<int>&);
    void SetAtomBondInfo();
    static void WarnBondLengthDefault(Atom::AtomicElementType,
                                      Atom::AtomicElementType,double);
    double GetBondLength(Atom::AtomicElementType, Atom::AtomicElementType);
    void GetBondsFromAtomCoords();
    void VisitAtom(int, int);
    int DetermineMolecules();
    void AtomDistance(int, int, int, std::set<int>&);
    void DetermineExcludedAtoms();
    void SetSolventInfo();

    void Mask_SelectDistance( Frame &, char*, bool, bool, double );
    void Mask_AND(char*,char*);
    void Mask_OR(char*,char*);
    void Mask_NEG(char*);
    void MaskSelectResidues(int, int, char *);
    // TODO: NameType Reference
    void MaskSelectResidues(NameType, char *);
    void MaskSelectElements( NameType, char* );
    void MaskSelectTypes( NameType, char* );
    void MaskSelectAtoms(int,int,char*);
    void MaskSelectAtoms(NameType,char*);
    bool ParseMask(Frame &, AtomMask &,bool);

    std::vector<int> SetupSequentialArray(std::vector<int>&, int, std::vector<int>&);
};
#endif
