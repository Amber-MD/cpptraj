#ifndef INC_TRAJ_PDBFILE_H
#define INC_TRAJ_PDBFILE_H
#include "TrajectoryIO.h"
#include "PDBtype.h"
// Class: Traj_PDBfile
/// TrajectoryIO class for reading coordinates from PDB files.
class Traj_PDBfile: public TrajectoryIO, PDBtype {
  public:
    /** PDBWRITEMODE: Indicate how the pdb should be written.
      *  SINGLE: Writing only a single frame.
      *  MODEL: Multiple frames written to the same file separated with 
      *         the MODEL keyword
      *  MULTI: Each frame written to a different file with name filename.frame
      */
    enum PDBWRITEMODE {SINGLE = 0, MODEL, MULTI};

    Traj_PDBfile();
    static TrajectoryIO* Alloc() { return (TrajectoryIO*)new Traj_PDBfile(); }
    // Traj_PDBfile-specfic functions
    void SetWriteMode(PDBWRITEMODE);
    void SetDumpq();

  private:
    int pdbAtom_;
    PDBWRITEMODE pdbWriteMode_;
    bool dumpq_; ///< If true, print charges in Occupancy column
    bool dumpr_; ///< If true, print radii in B-factor column.
    Topology *pdbTop_;

    std::vector<char> chainID_;
    char chainchar_;

    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(std::string const&, Topology*, TrajInfo&);
    int setupTrajout(std::string const&, Topology*, int, TrajInfo const&,bool);
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    void info();
    int processWriteArgs(ArgList&);
    int readVelocity(int, double*) { return 1; }
    int readIndices(int,int*) { return 1; }
    int processReadArgs(ArgList&) { return 0; }
};
#endif
