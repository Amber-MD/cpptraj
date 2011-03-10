#ifndef INC_AMBERPARM_H
#define INC_AMBERPARM_H
#include "PtrajFile.h" 

// Default size for atom and residue names, 4 + NULL
// NOTE: Also defined in ptrajmask.h
#define NAMESIZE 5

class AmberParm {
  public:
    // Default type for atom names, res names, atom types
    typedef char NAME[NAMESIZE];
  private:
    // ENUMERATED TYPE for TOPOLOGY VALUES
    enum topValues {
    //0       1       2      3       4       5       6       7      8       9
      NATOM,  NTYPES, NBONH, MBONA,  NTHETH, MTHETA, NPHIH,  MPHIA, NHPARM, NPARM, 
      NNB,    NRES,   NBONA, NTHETA, NPHIA,  NUMBND, NUMANG, NPTRA, NATYP,  NPHB, 
      IFPERT, NBPER,  NGPER, NDPER,  MBPER,  MGPER,  MDPER,  IFBOX, NMXRS, IFCAP,
      NEXTRA
    };
    // Enumerated type for Fortran Format
    enum FortranFormat {
      UNKNOWN_FFORMAT, F10I8, F5E16_8, F20a4, F12I6, F3I8
    };
    // Enumerated type for Fortran data type
    enum FortranType {
      UNKNOWN_FTYPE, FINT, FDOUBLE, FCHAR
    };
    // Contain data for SA LCPO
    struct SurfInfoType {
      double vdwradii;
      double P1;
      double P2;
      double P3;
      double P4;
    };
    typedef struct SurfInfoType SurfInfo;

    int *values;         // POINTERS

    int debug;
    // For determining fortran format
    FortranFormat fFormat;
    FortranType fType;
    int numCols;
    int width;
    const char *FormatString;
    int BufferSize;

    void AssignLCPO(SurfInfo *, double, double, double, double, double);
    void SetFormat(char *,int);
    void *getFlagFileValues(const char*,int);
    char *DataToBuffer(char*,const char *, int *, double *, NAME *, int N);
    int ReadParmMol2();
    int ReadParmAmber();
    int ReadParmPDB();

  public:

    PtrajFile File;
    char *parmName;       // Separate from File.filename in case of stripped parm
    int pindex;           // The index of this parm in the parmfilelist
    int parmFrames;       // For output, # of frames that will be read with this parm
    int outFrame;         // Output, # frames that have been written using this parm
    // Amber Parmtop
    int NbondsWithH();    // NBONH
    int NbondsWithoutH(); // MBONA
    int *bondsh;          // IBH/JBH/ICBH(NBONH)
    int *bonds;           // IB/JB/ICB(NBONA) NOTE: Using MBONA
    NAME *names;         // IGRAPH(NATOM)
    NAME *resnames;      // LBRES(NRES)
    NAME *types;         // ISYMBL(NATOM)
    char *ResidueName(int);
    int *resnums;         // IPRES(NRES) 
    int natom;            // NATOM
    int nres;             // NRES
    int ifbox;            // IFBOX
    int finalSoluteRes;   // IPTRES
    int molecules;        // NSPM
    int firstSolvMol;     // NSPSOL
    int *atomsPerMol;     // NSP(NSPM)
    double *mass;         // AMASS(NATOM)
    double *charge;       // CHARGE(NATOM)
    double *Box;          // OLDBETA, BOX(1), BOX(2), BOX(3)
    // From Ptraj
    char *solventMask;         // T for atoms in the solvent
    int solventMolecules;      // number of solvent molecules
    int *solventMoleculeStart; // pointer into solventMask for first atom of each solvent
    int *solventMoleculeStop;  // pointer into solventMask for last atom of each solvent
    int solventAtoms;          // number of solvent atoms
    // For SA calc
    SurfInfo *SurfaceInfo;

    AmberParm(int);
    ~AmberParm();
    void ResName(char *, int);
    int OpenParm(char *);
    int SetSurfaceInfo();
    int SetSolventInfo();
    void AtomInfo(int);
    void Info(char *);
    int atomToResidue(int);
    int atomToMolecule(int);
    int atomToSolventMolecule(int);
    AmberParm *modifyStateByMask(int *, int);
    AmberParm *modifyStateByMap(int *);
    int WriteAmberParm(); 
};
#endif
