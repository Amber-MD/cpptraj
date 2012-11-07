#ifndef INC_TRAJ_MOL2FILE_H
#define INC_TRAJ_MOL2FILE_H
#include "TrajectoryIO.h"
#include "Mol2File.h"
// Class: Traj_Mol2File
/// TrajecttoryIO class for reading coordinates from Mol2 files.
class Traj_Mol2File : public TrajectoryIO {
  public:
    /// Indicate how the mol2 file should be written.
    /** - SINGLE: Writing only a single frame
      * - MOL: Multiple frames written to the same file separated with
      *        a @<TRIPOS>MOLECULE section.
      * - MULTI: Each frame written to a different file with name filename.frame
      */
    enum MOL2WRITEMODE { NONE = 0, SINGLE, MOL, MULTI };

    Traj_Mol2File();
    static TrajectoryIO* Alloc() { return (TrajectoryIO*)new Traj_Mol2File(); }
  private:
    MOL2WRITEMODE mol2WriteMode_;
    Topology* mol2Top_;
    bool hasCharges_;
    std::vector<int> trajBonds_;
    Mol2File file_; 

    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(std::string const&, Topology*);
    int setupTrajout(std::string const&, Topology*, int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    void Info();
    int processWriteArgs(ArgList&);
    int readVelocity(int, double*) { return 1; }
    int readIndices(int,int*) { return 1; }
    int processReadArgs(ArgList&) { return 0; }
};
#endif
