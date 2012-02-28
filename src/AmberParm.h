#ifndef INC_AMBERPARM_H
#define INC_AMBERPARM_H
#include <vector>
#include "BoxType.h" 
// Name.h has definition for NAME 
#include "Name.h"
#include "Bonds.h"
#include "AtomMask.h"
// Class: AmberParm 
/** Hold all data pertaining to a molecular system (# atoms, atom names, 
  * etc). Can be read in from Amber Topology, PDB, Mol2, or Charmm PSF 
  * files (see ParmFile.cpp). The following parameters of AmberParm must 
  * always be set for the mask parser to function correctly:
  * - Arrays: names, resnames, resnums
  * - Variables: natom, nres 
  * Other actions may require additional info, such as bonding or surface area
  * terms. Note that in accordance with the rest of Amber, cpptraj expects 
  * atom/residue names etc (type defined in Name.h) to be uniformly 4 characters 
  * long, left-aligned, and it is up to the ReadParmXXX function to ensure this 
  * is the case. For example:
  * - Correct:   [CA  ]
  * - Incorrect: [CA]
  * - Incorrect: [ CA ]
  */
class AmberParm {
  private:
    int debug;
    /// Contain data for SA LCPO
    struct SurfInfo {
      double vdwradii;
      double P1;
      double P2;
      double P3;
      double P4;
    };
    void AssignLCPO(SurfInfo*,double,double,double,double,double);
    int numSoluteAtoms;

    // Set up solvent info
    bool IsSolventResname(NAME);
    int SetSolventInfo();
    bool hasSolventInfo;       ///< True if solvent information present

    // Routines to fill in missing Bond/Molecule info
    double *parmCoords;
    int DetermineMolecules();
    int SetupExcludedAtoms();
    BondInfo bondInfo;  ///< Class that holds bond info in a different format than bond arrays

    // Interface to atom mask parser in PtrajMask
    int SetupAtomMask(AtomMask &, double *, bool);

    // Amber Parmtop
    int *numex;         ///< NUMEX(NATOM)
    int *atype_index;   ///< IAC(NATOM)
    int ntypes;         ///< NTYPES
    int *at_num;
    int *NB_index;      ///< ICO(NTYPES*NTYPES)
    double *LJ_A;       ///< CN1(NTYPES*(NTYPES+1)/2)
    double *LJ_B;       ///< CN2(NTYPES*(NTYPES+1)/2)
    int nnb;            ///< NNB
    int *excludedAtoms; ///< INB(NNB)
    char *radius_set;   ///< TYPE 
    double *gb_radii;   ///< RBORN(NATOM)
    double *gb_screen;  ///< FS(NATOM)

    double *bond_rk;    ///< RK(NUMBND)
    double *bond_req;   ///< REQ(NUMBND)
    double *angle_tk;   ///< TK(NUMANG)
    double *angle_teq;  ///< TEQ(NUMANG)
    double *dihedral_pk;///< PK(NPTRA)
    double *dihedral_pn;///< PN(NPTRA)
    double *dihedral_phase; ///< PHASE(NPTRA)
    double *scee_scale; ///< ONE_SCEE(NPTRA)
    double *scnb_scale; ///< ONE_SCNB(NPTRA)
    double *solty;      ///< SOLTY(NATYP)
    int *anglesh;       ///< ITH/JTH/KTH/ICTH(NTHETH)
    int *angles;        ///< IT/JT/KT/ICT(NTHETA)
    int *dihedralsh;    ///< IPH/JPH/KPH/LPH/ICPH(NPHIH)
    int *dihedrals;     ///< IP/JP/KP/LP/ICP(NPHIA)
    double *asol;       ///< ASOL(NPHB)
    double *bsol;       ///< BSOL(NPHB)
    double *hbcut;      ///< HBCUT(NPHB)
    NAME *itree;        ///< ITREE(NATOM)
    int *join_array;    ///< JOIN(NATOM)
    int *irotat;        ///< IROTAT(NATOM)

    int numbnd;         ///< NUMBND: number of bond types
    int numang;         ///< NUMANG: number of angle types
    int numdih;         ///< NPTRA: number of dihedral types
    int NanglesWithH;   ///< NTHETH: number of angles containing hydrogen
    int NanglesWithoutH;///< MTHETA/(NTHETA): number of angles not containing hydrogen
    int NdihedralsWithH;///< NPHIH: number of dihedrals containing hydrogen
    int NdihedralsWithoutH; ///< MPHIA/(NPHIA): number of dihedrals not containing hydrogen
    int natyp;          ///< NATYP: number of atom types in parameter file (SOLTY)
    int nphb;           ///< NPHB: number of distinct 10-12 hydrogen bond pair types

    int NbondsWithH;    ///< NBONH
    int NbondsWithoutH; ///< MBONA
    int *bondsh;        ///< IBH/JBH/ICBH(NBONH)
    int *bonds;         ///< IB/JB/ICB(NBONA) NOTE: Using MBONA
    NAME *names;        ///< IGRAPH(NATOM)
    NAME *resnames;     ///< LBRES(NRES)
    NAME *types;        ///< ISYMBL(NATOM)
    int *resnums;       ///< IPRES(NRES) 
    int nres;           ///< NRES
    int finalSoluteRes; ///< IPTRES
    int molecules;      ///< NSPM
    int firstSolvMol;   ///< NSPSOL
    int *atomsPerMol;   ///< NSP(NSPM)
    double *charge;     ///< CHARGE(NATOM)

  public:
    char *parmfileName;   ///< Parm filename (full path)
    char *parmName;       ///< Parm name, set to base filename on reads 
    int pindex;           ///< The index of this parm in the parmfilelist
    int parmFrames;       ///< For output, # of frames that will be read with this parm

    int natom;          ///< NATOM
    double *mass;       ///< AMASS(NATOM)

    // From Ptraj
    // NOTE: Eventually want to make this private, or eliminate it entirely.
    //       Kept for now since some actions in ptraj_actions.c require it.
    char *solventMask;         ///< T for atoms in the solvent
    int solventMolecules;      ///< number of solvent molecules
    int *solventMoleculeStart; ///< pointer into solventMask for first atom of each solvent
    int *solventMoleculeStop;  ///< pointer into solventMask for last atom of each solvent
    int solventAtoms;          ///< number of solvent atoms

    int Nres()           { return nres;           }
    int FirstSolventMol(){ return firstSolvMol;   }
    int Nmol()           { return molecules;      }
    bool HasSolventInfo(){ return hasSolventInfo; }

    NAME *AtomNames_ptr()    { return names;       } ///< Returns array IGRAPH
    NAME *ResidueNames_ptr() { return resnames;    } ///< Returns array LBRES
    NAME *AtomTypes_ptr()    { return types;       } ///< Returns array ISYMBL
    int *AtomsPerMol_ptr()   { return atomsPerMol; } ///< Returns array NSP
    int *ResAtomNums_ptr()   { return resnums;     } ///< Returns array IPRES
    double *Charges_ptr()    { return charge;      } ///< Returns array CHARGE
    double *GB_radii_ptr()   { return gb_radii;    } ///< Returns array RBORN

    int SetupIntegerMask(AtomMask &, double *);
    int SetupCharMask(AtomMask &, double *);
    int NumMoleculesInMask(AtomMask &);

    int SetupExcludedAtomsList(AtomMask &, std::vector< std::vector<int> > &);
    int GetLJparam(double *, double *, int, int);
    int GetBondParamIdx(int, double *, double *);
    double GetBondedCutoff(int, int);
    //int GetBondParam(double *, double *, int, int);
    int SetCharges(double*);
    /// Set the input array with atom charges in Amber format (pre-multiplied by 18.2223)
    int AmberChargeArray(std::vector<double>&);
    double AtomCharge(int);
    int AtomsPerMol(int );

    double Box[6];      ///< X, Y, Z, alpha, beta, gamma 
    BoxType boxType;    ///< None, Orthogonal, Non-orthogonal

    // For SA calc
    SurfInfo *SurfaceInfo;

    AmberParm();
    ~AmberParm();

    void SetDebug(int);

    void ResName(char *, int);
    void ResAtomName(char*, int);
    char *ResidueName(int);
    int FindAtomInResidue(int, char *);
    int ResAtomRange(int, int *, int *);
    int FinalSoluteRes();

    char *AtomName(int);
    bool AtomNameIs(int, const char *);
    bool AtomElementIs(int, AtomicElementType);

    int SetSurfaceInfo();

    int CommonSetup(bool, bool);
    void SetParmName(char *, char *);

    void AtomInfo(int);
    void ParmInfo(std::string&);
    void Summary();
    void PrintBondInfo();
    void PrintMoleculeInfo();
    void PrintResidueInfo();

    int atomToResidue(int);
    int atomToMolecule(int);
    int atomToSolventMolecule(int);

    void GetBondsFromCoords(double*);
    int BondArray(std::vector<int> &);
    int BondArrayWithParmIdx(std::vector<int> &);
    int SetupBondInfo();
    int GetBondedAtomIdx(int, const char *);
    int GetBondedHatoms(int,std::vector<int>&);
    int MaskOfAtomsAroundBond(int, int, std::vector<char>&);

    void ResetBondInfo(); 
    int AddBond(int, int, int);

    AmberParm *modifyStateByMask(std::vector<int>&, char*);
    AmberParm *modifyStateByMap(int *);

    // Friend classes for reading/writing parm info
    friend class AmberParmFile;
    friend class PdbParmFile;
    friend class Mol2ParmFile;
    friend class CharmmPsfParmFile; 
};
#endif
