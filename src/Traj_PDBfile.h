#ifndef INC_TRAJ_PDBFILE_H
#define INC_TRAJ_PDBFILE_H
#include "TrajectoryIO.h"
// Class: PDBfile
/// TrajectoryIO class for reading coordinates from PDB files.
class PDBfile: public TrajectoryIO {
  public:
    /** PDBWRITEMODE: Indicate how the pdb should be written.
      *  SINGLE: Writing only a single frame.
      *  MODEL: Multiple frames written to the same file separated with 
      *         the MODEL keyword
      *  MULTI: Each frame written to a different file with name filename.frame
      */
    enum PDBWRITEMODE {SINGLE = 0, MODEL, MULTI};

    PDBfile();
    // PDBfile-specfic functions
    void SetWriteMode(PDBWRITEMODE);
    void SetDumpq();

  private:
    static const size_t BUF_SIZE;

    int pdbAtom_;
    PDBWRITEMODE pdbWriteMode_;
    bool dumpq_; ///< If true, print charges in Occupancy column
    bool dumpr_; ///< If true, print radii in B-factor column.
    // The following are only required for writes and are set in setupTrajout 
    NAME *pdbAtomNames_; 
    NAME *trajResNames_;
    int *trajAtomsPerMol_;
    int *trajResNums_;
    double *trajCharges_;
    double *trajRadii_;

    std::vector<char> chainID_;
    char chainchar_;

    // Inherited functions
    int setupTrajin(AmberParm*);
    int setupTrajout(AmberParm*);
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    void info();
    int processWriteArgs(ArgList*);
};
#endif
