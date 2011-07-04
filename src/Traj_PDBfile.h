#ifndef INC_TRAJ_PDBFILE_H
#define INC_TRAJ_PDBFILE_H
#include "TrajectoryIO.h"
/// Class: PDBfile
/// TrajectoryIO class for reading coordinates from PDB files.
class PDBfile: public TrajectoryIO {
  public:
    // PDBWRITEMODE: Indicate how the pdb should be written.
    //   SINGLE: Writing only a single frame.
    //   MODEL: Multiple frames written to the same file separated with 
    //          the MODEL keyword
    //   MULTI: Each frame written to a different file with name filename.frame
    enum PDBWRITEMODE {SINGLE = 0, MODEL, MULTI};
  private:
    char buffer[256];
    int pdbAtom;
    PDBWRITEMODE pdbWriteMode;
    bool dumpq; // If true, print charges in Occupancy column
    NAME *pdbAtomNames; // Read in during setupRead for checking against
                        // parm names via CheckPdbNames, or set for writes
                        // using SetParmInfo
    // The following are only required for writes and are set in SetParmInfo
    NAME *trajResNames;
    int *trajAtomsPerMol;
    int *trajResNums;
    double *trajCharges;
    double *trajRadii;

    // Inherited functions
    int setupRead(int);
    int setupWrite(int);
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    void info();

  public:
    PDBfile();
    ~PDBfile();
    // PDBfile-specfic functions
    void SetWriteMode(PDBWRITEMODE);
    void SetDumpq() { dumpq = true; }
    int CheckPdbNames(NAME*);
    void NumFramesToWrite(int);
    void SetParmInfo(NAME*, NAME*, int *, int *, double *, double *);
};
#endif
