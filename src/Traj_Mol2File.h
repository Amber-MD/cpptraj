#ifndef INC_TRAJ_MOL2FILE_H
#define INC_TRAJ_MOL2FILE_H
#include "TrajectoryIO.h"
/// Class: Mol2File
/// TrajecttoryIO class for reading coordinates from Mol2 files.
class Mol2File : public TrajectoryIO {
  public:
    // MOL2WRITEMODE: Indicate how the mol2 file should be written.
    //   SINGLE: Writing only a single frame
    //   MOL: Multiple frames written to the same file separated with
    //        a @<TRIPOS>MOLECULE section.
    //   MULTI: Each frame written to a different file with name filename.frame
    enum MOL2WRITEMODE { SINGLE = 0, MOL, MULTI };
  private:
    int mol2atom;
    int mol2bonds;
    MOL2WRITEMODE mol2WriteMode;

    // The following are only required for writes and are set in setupWrite 
    int trajnres;
    int trajnbondsh;
    int trajnbonds;
    AmberParm::NAME *trajAtomNames; 
    AmberParm::NAME *trajTypes;
    AmberParm::NAME *trajResNames; 
    int *trajResNums;
    double *trajCharges;
    int *trajBonds;
    int *trajBondsh;

    // Inherited functions
    int setupRead(AmberParm *);
    int setupWrite(AmberParm *);
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    void info();
    int processWriteArgs(ArgList *);
  public :
    Mol2File();
    ~Mol2File();
    // Mol2-specific functions
    void SetWriteMode(MOL2WRITEMODE);
};
#endif
