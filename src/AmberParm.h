#ifndef INC_AMBERPARM_H
#define INC_AMBERPARM_H
#include "CpptrajFile.h"
#include "BoxType.h" 
// Name.h has definition for NAME 
#include "Name.h"
// Class: AmberParm 
/** Hold all data pertaining to a molecular system (# atoms, atom names, 
  * etc). Can be read in from Amber Topology, PDB, or Mol2 files (implemented 
  * in the ReadParmXXX functions). The following parameters of AmberParm must 
  * always be set:
  *   1. Variables: natom, nres, boxType 
  *   2.    Arrays: names, resnames, resnums
  */
class AmberParm {
  private:
    int debug;
    /// ENUMERATED TYPE for AMBER TOPOLOGY VALUES
    enum topValues {
    //0       1       2      3       4       5       6       7      8       9
      NATOM,  NTYPES, NBONH, MBONA,  NTHETH, MTHETA, NPHIH,  MPHIA, NHPARM, NPARM, 
      NNB,    NRES,   NBONA, NTHETA, NPHIA,  NUMBND, NUMANG, NPTRA, NATYP,  NPHB, 
      IFPERT, NBPER,  NGPER, NDPER,  MBPER,  MGPER,  MDPER,  IFBOX, NMXRS, IFCAP,
      NEXTRA
    };

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

    // Parm format readers
    int ReadParmMol2(CpptrajFile *);
    void SetParmFromValues(int*,bool);
    int ReadParmOldAmber(CpptrajFile *);
    int ReadParmAmber(CpptrajFile *);
    int SetAtomsPerMolPDB(int);
    int ReadParmPDB(CpptrajFile *);
    int ReadParmPSF(CpptrajFile *);

    // Routines to fill in missing Bond/Molecule info
    void GetBondsFromCoords();
    double *parmCoords;
    int DetermineMolecules();
    int SetupExcludedAtoms();

    int *numex;         ///< NUMEX(NATOM)
    int *atype_index;   ///< IAC(NATOM)
    int ntypes;         ///< NTYPES
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

  public:
    char *parmfileName;   ///< Parm filename (full path)
    char *parmName;       ///< Parm name, set to base filename on reads 
    int pindex;           ///< The index of this parm in the parmfilelist
    int parmFrames;       ///< For output, # of frames that will be read with this parm

    // Amber Parmtop
    int NbondsWithH;    ///< NBONH
    int NbondsWithoutH; ///< MBONA
    int *bondsh;        ///< IBH/JBH/ICBH(NBONH)
    int *bonds;         ///< IB/JB/ICB(NBONA) NOTE: Using MBONA
    NAME *names;        ///< IGRAPH(NATOM)
    NAME *resnames;     ///< LBRES(NRES)
    NAME *types;        ///< ISYMBL(NATOM)
    int *resnums;       ///< IPRES(NRES) 
    int natom;          ///< NATOM
    int nres;           ///< NRES
    int finalSoluteRes; ///< IPTRES
    int molecules;      ///< NSPM
    int firstSolvMol;   ///< NSPSOL
    int *atomsPerMol;   ///< NSP(NSPM)
    double *mass;       ///< AMASS(NATOM)
    double *charge;     ///< CHARGE(NATOM)

    double *GB_radii() { return gb_radii; } ///< RBORN(NATOM)

    int NumExcludedAtoms(int);
    int Natex(int);
    int GetLJparam(double *, double *, int, int);
    int GetBondParamIdx(int, double *, double *);
    //int GetBondParam(double *, double *, int, int);
    int SetCharges(double*);

    double Box[6];      ///< X, Y, Z, alpha, beta, gamma 
    BoxType boxType;    ///< None, Orthogonal, Non-orthogonal

    // From Ptraj
    char *solventMask;         ///< T for atoms in the solvent
    int solventMolecules;      ///< number of solvent molecules
    int *solventMoleculeStart; ///< pointer into solventMask for first atom of each solvent
    int *solventMoleculeStop;  ///< pointer into solventMask for last atom of each solvent
    int solventAtoms;          ///< number of solvent atoms

    // For SA calc
    SurfInfo *SurfaceInfo;

    AmberParm();
    ~AmberParm();

    void SetDebug(int);

    void ResName(char *, int);
    void ResAtomName(char*, int);
    char *ResidueName(int);
    int FindAtomInResidue(int, char *);

    int SetSurfaceInfo();

    int OpenParm(char *,bool,bool);

    void AtomInfo(int);
    void ParmInfo();
    void Summary();
    void PrintBondInfo();
    void PrintMoleculeInfo();

    int atomToResidue(int);
    int atomToMolecule(int);
    int atomToSolventMolecule(int);

    void ResetBondInfo(); 
    int AddBond(int, int, int);

    AmberParm *modifyStateByMask(int *, int, char*);
    AmberParm *modifyStateByMap(int *);

    int WriteAmberParm(char*); 
};
#endif
