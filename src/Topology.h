#ifndef INC_TOPOLOGY_H
#define INC_TOPOLOGY_H
#include <string>
#include "Atom.h"
#include "Residue.h"
#include "Molecule.h"
#include "AtomMask.h"
#include "Box.h"
#include "CoordFrame.h"
#include "Frame.h"
#include "ParameterTypes.h"
// Class: Topology
/// Hold information for all atoms
class Topology {
  public:
    Topology();
    ~Topology();
    // ----- Set internal variables -----
    void SetDebug(int);
    void SetHasCoordinates(); // TODO: Replace with pass-in of CoordFrame
    void SetParmName(const char*);
    void SetParmName(std::string&);
    void SetPindex(int);
    void SetReferenceCoords( Frame* ); // TODO: Pass in frame reference
    // ----- Return internal variables -----
    int FinalSoluteRes();
    // ---- Atom-specific routines -----
    typedef std::vector<Atom>::const_iterator atom_iterator;
    inline atom_iterator begin() const {
      return atoms_.begin();
    }
    inline atom_iterator end() const {
      return atoms_.end();
    }
    atom_iterator ResAtomStart(int) const;
    atom_iterator ResAtomEnd(int) const;
    atom_iterator MolAtomStart(int) const;
    atom_iterator MolAtomEnd(int) const;
    const Atom &operator[](int); // NOTE: Inline?
    // ----- Residue-specific routines -----
    typedef std::vector<Residue>::const_iterator res_iterator;
    inline res_iterator ResStart() const {
      return residues_.begin();
    }
    inline res_iterator ResEnd() const {
      return residues_.end();
    }
    const Residue &Res(int); 
    // ----- Molecule-specific routines -----
    typedef std::vector<Molecule>::const_iterator mol_iterator;
    inline mol_iterator MolStart() const {
      return molecules_.begin();
    }
    inline mol_iterator MolEnd() const {
      return molecules_.end();
    }
    mol_iterator SolventStart() const;
    mol_iterator SolventEnd() const;
    // ----- Bond-specific routines -----
    typedef std::vector<int>::const_iterator bond_iterator;
    inline bond_iterator BondsStart() const {
      return bonds_.begin();
    }
    inline bond_iterator BondsEnd() const {
      return bonds_.end();
    }
    inline bond_iterator BondsH_Start() const {
      return bondsh_.begin();
    }
    inline bond_iterator BondsH_End() const {
      return bondsh_.end();
    }
    // ----- Misc routines -----
    int ResAtomRange(int, int *, int *);
    char *ResidueName(int); // TODO: Make obsolete
    std::string ResAtomName(int);
    int FindAtomInResidue(int, NameType);
    int FindResidueMaxNatom();
    int SoluteAtoms();
    double *Mass();
    // ----- Print topology info -----
    void Summary();
    void ParmInfo();
    void PrintAtomInfo(const char*);
    void PrintBondInfo();
    // ----- Return Internal Variables -----
    inline int Pindex() {
      return pindex_;
    }
    inline int Natom() {
      return (int)atoms_.size();
    }
    inline int Nres() {
      return (int)residues_.size();
    }
    inline int Nmol() {
      return (int)molecules_.size();
    }
    inline int FirstSolventMol() {
      return firstSolventMol_;
    }
    inline int Nsolvent() {
      return NsolventMolecules_;
    }
    inline int Nframes() {
      return nframes_;
    }
    inline void IncreaseFrames(int frames) {
      nframes_ += frames;
    }
    inline const char *ResName(int resnum) {
      return residues_[resnum].c_str();
    }
    inline int Mol_FirstRes(int mol) {
      return molecules_[mol].FirstRes();
    }
    inline const char *c_str() {
      return parmName_.c_str();
    }
    inline std::string &ParmName() {
      return parmName_;
    }
    // ----- Routines to Access/Modify Box info -----
    inline Box::BoxType BoxType() {
      return box_.Type();
    }
    inline bool BoxIsTruncOct() { // NOTE: Should this not be inlined?
      return (box_.AmberIfbox() == 2);
    }
    inline void SetNoBox() {
      box_.SetNoBox();
    }
    inline void SetBoxAngles(double *abgIn) {
      box_.SetAngles( abgIn );
    }
    inline void BoxCoords(double *boxOut) {
      box_.ToDouble(boxOut);
    }
    // ----- PDB/Mol2 etc setup routines -----
    void AddAtom(Atom, Residue);
    void StartNewMol();
    // ----- Amber setup routines -----
    int CreateAtomArray(std::vector<NameType>&, std::vector<double>&,
                        std::vector<int>&, std::vector<double>&,
                        std::vector<int>&, std::vector<NameType>&,
                        std::vector<double>&,std::vector<double>&, 
                        std::vector<NameType>&, std::vector<int>&);
    int CreateMoleculeArray(std::vector<int> &,Box,int,int);
    int SetBondInfo(std::vector<int> &, std::vector<int> &,
                    std::vector<double>&,std::vector<double>&);

    void CommonSetup(bool,bool);

    void ClearBondInfo();
    void AddBond(int,int);
    // ----- Mask Routines -----
    inline bool SetupIntegerMask(AtomMask &mask) {
      return ParseMask(refCoords_, mask, true);
    }
    inline bool SetupCharMask(AtomMask &mask) {
      return ParseMask(refCoords_, mask, false);
    }
    bool SetupIntegerMask(AtomMask &, Frame &);
    bool SetupCharMask(AtomMask &, Frame &);

    Topology *modifyStateByMask(AtomMask &, const char *);

  private:
    std::vector<Atom> atoms_;
    std::vector<Residue> residues_;
    std::vector<Molecule> molecules_;
    std::string parmName_;

    std::vector<int> bonds_;
    std::vector<int> bondsh_;
    std::vector<ParmBondType> BondParm;
    std::vector<int> excludedAtoms_;
    std::vector<int> numex_;

    Box box_;
    CoordFrame refCoords_;

    int debug_;
    bool hasCoordinates_;
    int topology_error_;
    int firstSolventMol_;
    int NsolventMolecules_;
    int finalSoluteRes_;
    int pindex_;
    int nframes_;
    double *massptr_;

    void SetAtomBondInfo();
    double GetBondedCut(Atom::AtomicElementType, Atom::AtomicElementType);
    void GetBondsFromAtomCoords();
    void VisitAtom(int, int);
    void DetermineMolecules();
    void AtomDistance(int, int);
    void DetermineExcludedAtoms();
    void SetSolventInfo();

    void Mask_SelectDistance( CoordFrame &, char*, bool, bool, double );
    void Mask_AND(char*,char*);
    void Mask_OR(char*,char*);
    void Mask_NEG(char*);
    void MaskSelectResidues(int, int, char *);
    // TODO: NameType Reference
    void MaskSelectResidues(NameType, char *);
    void MaskSelectAtoms(int,int,char*);
    void MaskSelectAtoms(NameType,char*);
    bool ParseMask(CoordFrame &, AtomMask &,bool);

    std::vector<int> SetupSequentialArray(std::vector<int>&, int, std::vector<int>&);


};
#endif
