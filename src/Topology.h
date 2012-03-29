#ifndef INC_TOPOLOGY_H
#define INC_TOPOLOGY_H
#include <string>
#include "Atom.h"
#include "Residue.h"
#include "Molecule.h"
#include "AtomMask.h"
#include "Box.h"
// Class: Topology
/// Hold information for all atoms
class Topology {
  public:
    Topology();

    void SetDebug(int);
    void SetHasCoordinates();
    void SetParmName(const char*);
    void SetParmName(std::string&);
    void SetPindex(int);

    typedef std::vector<Atom>::const_iterator atom_iterator;
    atom_iterator begin() const;
    atom_iterator end() const;
    atom_iterator ResStart(int) const;
    atom_iterator ResEnd(int) const;
    atom_iterator MolStart(int) const;
    atom_iterator MolEnd(int) const;

    //typedef std::vector<Residue>::const_iterator res_iterator;
    //res_iterator beginRes() const;
    //res_iterator endRes() const;

    void Summary();
    void ParmInfo();
    void PrintAtomInfo(const char*);
    void PrintBondInfo();

    inline int Pindex() {
      return pindex_;
    }
    inline int Nres() {
      return (int)residues_.size();
    }
    inline int Nmol() {
      return (int)molecules_.size();
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

    void AddAtom(Atom, Residue);
    void StartNewMol();

    int CreateAtomArray(std::vector<NameType>&, std::vector<double>&,
                        std::vector<int>&, std::vector<double>&,
                        std::vector<int>&,
                        std::vector<double>&,std::vector<double>&, 
                        std::vector<NameType>&, std::vector<int>&);
    int CreateMoleculeArray(std::vector<int> &,Box,int,int);
    int SetBondInfo(std::vector<int> &, std::vector<int> &);

    void CommonSetup(bool,bool);

    int ParseMask(AtomMask &,bool);
  private:
    static const char AtomicElementName[][3];

    std::vector<Atom> atoms_;
    std::vector<Residue> residues_;
    std::vector<Molecule> molecules_;
    std::string parmName_;

    std::vector<int> bonds_;
    std::vector<int> bondsh_;
    std::vector<int> excludedAtoms_;
    std::vector<int> numex_;

    Box box;

    int debug_;
    bool hasCoordinates_;
    int topology_error_;
    int firstSolventMol_;
    int finalSoluteRes_;
    int pindex_;

    double GetBondedCut(Atom::AtomicElementType, Atom::AtomicElementType);
    void GetBondsFromAtomCoords();
    void AddBond(int,int);
    void VisitAtom(int, int);
    void DetermineMolecules();
    void AtomDistance(int, int);
    void DetermineExcludedAtoms();
    void SetSolventInfo();

    void Mask_SelectDistance( char*, bool, bool, double );
    void Mask_AND(char*,char*);
    void Mask_OR(char*,char*);
    void Mask_NEG(char*);
    void MaskSelectResidues(int, int, char *);
    // TODO: NameType Reference
    void MaskSelectResidues(NameType, char *);
    void MaskSelectAtoms(int,int,char*);
    void MaskSelectAtoms(NameType,char*);

    //std::vector<size_t> resnums_;
    //std::vector<NameType> resnames_;
    //int nres_; ///< Actual # of residues, 1 less than size of resnums
/*    /// Hold information for a bond
    struct BondType {
      double rk;
      double req;
      //double scale; In ptraj - necessary?
      int at1;
      int at2;
      int bidx;
    };
    /// NBOND: Hold information for bonds without hydrogen
    std::vector<BondType> Bond;
    /// NBONDH: Hold information for bonds with hydrogen
    std::vector<BondType> BondH;

    /// Hold information for an angle
    struct AngleType {
      double tk;
      double teq;
      int at1;
      int at2;
      int at3;
      int aidx;
    };
    /// NTHETA: Hold information for angles without hydrogen
    std::vector<AngleType> Angle;
    /// NTHETH: Hold information for angles with hydrogen
    std::vector<AngleType> AngleH;

    /// Hold information for a dihedral
    struct DihedralType {
      double pk;
      double pn;
      double phase;
      int at1;
      int at2;
      int at3;
      int at4;
      int didx;
    };

    /// Hold information for LCPO Surface Area calc
    struct LcpoType {
      double vdwradii;
      double P1;
      double P2;
      double P3;
      double P4;
    };
    /// NATOM: Hold LCPO info for all atoms
    std::vector<LcpoType> Lcpo;
*/

};
#endif
