#ifndef INC_AMBERPARM_H
#define INC_AMBERPARM_H
/// Class: AmberParm 
/// Hold all data pertaining to a molecular system (# atoms, atom names, 
/// etc). Can be read in from Amber Topology, PDB, or Mol2 files (implemented 
/// in the ReadParmXXX functions). The following parameters of AmberParm must 
/// always be set:
///   1. Variables: natom, nres, and ifbox
///   2.    Arrays: names, resnames, resnums
#include "PtrajFile.h"
#include "BoxType.h" 
// NAMESIZE: Default size for atom and residue names, 5 + NULL.
// Amber atom/residue names are 4, but some mol2 atom types are larger.
#include "Name.h"
class AmberParm {
  private:
    int debug;
    // ENUMERATED TYPE for AMBER TOPOLOGY VALUES
    enum topValues {
    //0       1       2      3       4       5       6       7      8       9
      NATOM,  NTYPES, NBONH, MBONA,  NTHETH, MTHETA, NPHIH,  MPHIA, NHPARM, NPARM, 
      NNB,    NRES,   NBONA, NTHETA, NPHIA,  NUMBND, NUMANG, NPTRA, NATYP,  NPHB, 
      IFPERT, NBPER,  NGPER, NDPER,  MBPER,  MGPER,  MDPER,  IFBOX, NMXRS, IFCAP,
      NEXTRA
    };

    // Contain data for SA LCPO
    struct SurfInfo {
      double vdwradii;
      double P1;
      double P2;
      double P3;
      double P4;
    };
    void AssignLCPO(SurfInfo*,double,double,double,double,double);

    // Set up solvent info
    bool IsSolventResname(NAME);
    int SetSolventInfo();

    // Parm format readers
    int ReadParmMol2(PtrajFile *);
    int ReadParmAmber(PtrajFile *);
    int SetAtomsPerMolPDB(int);
    int ReadParmPDB(PtrajFile *);

  public:
    char *parmfileName;   // Parm filename (full path)
    char *parmName;       // Parm name, set to base filename on reads 
    int pindex;           // The index of this parm in the parmfilelist
    int parmFrames;       // For output, # of frames that will be read with this parm
    //int outFrame;         // Output, # frames that have been written using this parm

    // Amber Parmtop
    int NbondsWithH;    // NBONH
    int NbondsWithoutH; // MBONA
    int *bondsh;        // IBH/JBH/ICBH(NBONH)
    int *bonds;         // IB/JB/ICB(NBONA) NOTE: Using MBONA
    NAME *names;        // IGRAPH(NATOM)
    NAME *resnames;     // LBRES(NRES)
    NAME *types;        // ISYMBL(NATOM)
    int *resnums;       // IPRES(NRES) 
    int natom;          // NATOM
    int nres;           // NRES
    int finalSoluteRes; // IPTRES
    int molecules;      // NSPM
    int firstSolvMol;   // NSPSOL
    int *atomsPerMol;   // NSP(NSPM)
    double *mass;       // AMASS(NATOM)
    double *charge;     // CHARGE(NATOM)
    double Box[6];      // X, Y, Z, alpha, beta, gamma 
    BoxType boxType;    // None, Orthogonal, Non-orthogonal

    // From Ptraj
    char *solventMask;         // T for atoms in the solvent
    int solventMolecules;      // number of solvent molecules
    int *solventMoleculeStart; // pointer into solventMask for first atom of each solvent
    int *solventMoleculeStop;  // pointer into solventMask for last atom of each solvent
    int solventAtoms;          // number of solvent atoms

    // For SA calc
    SurfInfo *SurfaceInfo;

    AmberParm();
    ~AmberParm();

    void SetDebug(int);

    void ResName(char *, int);
    void ResAtomName(char*, int);
    char *ResidueName(int);

    int SetSurfaceInfo();

    int OpenParm(char *);

    void AtomInfo(int);
    void ParmInfo();
    void Summary();

    int atomToResidue(int);
    int atomToMolecule(int);
    int atomToSolventMolecule(int);

    void ResetBondInfo(); 
    int AddBond(int, int, int);

    AmberParm *modifyStateByMask(int *, int);
    AmberParm *modifyStateByMap(int *);

    int WriteAmberParm(char*); 
};
#endif
