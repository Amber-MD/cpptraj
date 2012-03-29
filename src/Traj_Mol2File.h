#ifndef INC_TRAJ_MOL2FILE_H
#define INC_TRAJ_MOL2FILE_H
#include "TrajectoryIO.h"
#include "Mol2File.h"
// Class: Mol2File
/// TrajecttoryIO class for reading coordinates from Mol2 files.
class Mol2File : public TrajectoryIO, Mol2File {
  public:
    /// Indicate how the mol2 file should be written.
    /** - SINGLE: Writing only a single frame
      * - MOL: Multiple frames written to the same file separated with
      *        a @<TRIPOS>MOLECULE section.
      * - MULTI: Each frame written to a different file with name filename.frame
      */
    enum MOL2WRITEMODE { SINGLE = 0, MOL, MULTI };

    Mol2File();
    // Mol2-specific functions
    void SetWriteMode(MOL2WRITEMODE);
  private:
    int mol2atom_;
    int mol2bonds_;
    MOL2WRITEMODE mol2WriteMode_;

    // The following are only required for writes and are set in setupTrajout 
    int trajnres_;
    NAME *trajAtomNames_;
    NAME *trajTypes_;
    NAME *trajResNames_; 
    int *trajResNums_;
    double *trajCharges_;
    std::vector<int> trajBonds_;

    // Inherited functions
    int setupTrajin(Topology *);
    int setupTrajout(Topology *);
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    void info();
    int processWriteArgs(ArgList *);
};
#endif
