#ifndef INC_TRAJ_MOL2FILE_H
#define INC_TRAJ_MOL2FILE_H
#include "TrajectoryIO.h"
#include "Mol2File.h"
// Class: Traj_Mol2File
/// TrajecttoryIO class for reading coordinates from Mol2 files.
class Traj_Mol2File : public TrajectoryIO, Mol2File {
  public:
    /// Indicate how the mol2 file should be written.
    /** - SINGLE: Writing only a single frame
      * - MOL: Multiple frames written to the same file separated with
      *        a @<TRIPOS>MOLECULE section.
      * - MULTI: Each frame written to a different file with name filename.frame
      */
    enum MOL2WRITEMODE { SINGLE = 0, MOL, MULTI };

    Traj_Mol2File();
    // Mol2-specific functions
    void SetWriteMode(MOL2WRITEMODE);
  private:
    MOL2WRITEMODE mol2WriteMode_;
    Topology *mol2Top_;
    bool hasCharges_;
    std::vector<int> trajBonds_; 

    // Inherited functions
    bool ID_TrajFormat();
    int setupTrajin(Topology *);
    int setupTrajout(Topology *,int);
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    void info();
    int processWriteArgs(ArgList *);
};
#endif
