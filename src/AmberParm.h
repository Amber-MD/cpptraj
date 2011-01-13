#ifndef INC_AMBERPARM_H
#define INC_AMBERPARM_H
#include "PtrajFile.h" 

class AmberParm {
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
  void AssignLCPO(SurfInfo *, double, double, double, double, double);

  int *values;
  int *bonds;
  int debug;
  // For determining fortran format
  FortranFormat fFormat;
  FortranType fType;
  int numCols;
  int width;
  const char *FormatString;
  int BufferSize;

  void SetFormat(char *,int);
  void *getFlagFileValues(const char*,int);
  char *DataToBuffer(char*,const char *, int *, double *, char **, int N);

  public:
  int Nbonh() { return values[NBONH]; }
  int *bondsh;
  PtrajFile File;
  char *parmName;     // Separate from File.filename in case of stripped parm
  char **names;
  char **types;
  char **resnames;
  int *resnums;       // IPRES 
  int natom;
  int nres;
  int ifbox;
  int finalSoluteRes; // IPTRES
  int molecules;      // NSPM
  int firstSolvMol;   // NSPSOL
  int *atomsPerMol;   // NSP
  double *mass;
  double *charge;
  double *Box;
  // From Ptraj
  char *solventMask;         // T for atoms in the solvent
  int solventMolecules;      // number of solvent molecules
  int *solventMoleculeStart; // pointer into solventMask for first atom of each solvent
  int *solventMoleculeStop;  // pointer into solventMask for last atom of each solvent
  int solventAtoms;          // number of solvent atoms

  int pindex;       // The index of this parm in the parmfilelist
  int parmFrames;   // For output, # of frames that will be read with this parm
  int outFrame;     // Output, # frames that have been written using this parm
  // For SA calc
  SurfInfo *SurfaceInfo;

  AmberParm(int);
  ~AmberParm();
  void ResName(char *, int);
  int OpenParm(char *);
  int SetSurfaceInfo();
  int SetSolventInfo();
  int ReadParmAmber();
  int ReadParmPDB();
  void ParmInfo(char *);
  void Info(char *);
  char *mask(char *);
  char *mask(char *, double *);
  int atomToResidue(int);
  int atomToMolecule(int);
  int atomToSolventMolecule(int);
  AmberParm *modifyStateByMask(int *, int);
  int WriteAmberParm(); 
};
#endif
