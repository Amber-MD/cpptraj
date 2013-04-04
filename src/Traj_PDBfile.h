#ifndef INC_TRAJ_PDBFILE_H
#define INC_TRAJ_PDBFILE_H
#include "TrajectoryIO.h"
#include "PDBfile.h"
// Class: Traj_PDBfile
/// TrajectoryIO class for reading coordinates from PDB files.
class Traj_PDBfile: public TrajectoryIO {
  public:
    /** PDBWRITEMODE: Indicate how the pdb should be written.
      *  SINGLE: Writing only a single frame.
      *  MODEL: Multiple frames written to the same file separated with 
      *         the MODEL keyword
      *  MULTI: Each frame written to a different file with name filename.frame
      */
    enum PDBWRITEMODE {NONE = 0, SINGLE, MODEL, MULTI};

    Traj_PDBfile();
    static TrajectoryIO* Alloc() { return (TrajectoryIO*)new Traj_PDBfile(); }
  private:
    int pdbAtom_;
    int ter_num_; ///< Amount to increment atom number for TER
    PDBWRITEMODE pdbWriteMode_;
    bool dumpq_; ///< If true, print charges in Occupancy column
    bool dumpr_; ///< If true, print radii in B-factor column.
    Topology *pdbTop_;
    PDBfile file_;

    std::vector<char> chainID_;
    char chainchar_;

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
