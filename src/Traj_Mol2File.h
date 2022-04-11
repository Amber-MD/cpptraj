#ifndef INC_TRAJ_MOL2FILE_H
#define INC_TRAJ_MOL2FILE_H
#include "TrajectoryIO.h"
#include "Mol2File.h"
/// TrajectoryIO class for reading coordinates from Mol2 files.
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
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_Mol2File(); }
    static void WriteHelp();
  private:
    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(FileName const&, Topology*);
    int setupTrajout(FileName const&, Topology*, CoordinateInfo const&,int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,Frame&);
    int writeFrame(int,Frame const&);
    void Info();
    int processWriteArgs(ArgList&, DataSetList const&);
    int readVelocity(int, Frame&) { return 1; }
    int readForce(int, Frame&)    { return 1; }
    int processReadArgs(ArgList&) { return 0; }
#   ifdef MPI
    // Parallel functions
    int parallelOpenTrajout(Parallel::Comm const&);
    int parallelSetupTrajout(FileName const&, Topology*, CoordinateInfo const&,
                             int, bool, Parallel::Comm const&);
    int parallelWriteFrame(int, Frame const&);
    void parallelCloseTraj() {}
#   endif
    MOL2WRITEMODE mol2WriteMode_;
    Topology* mol2Top_;
    std::string ac_filename_;
    std::string bc_filename_;
    int currentSet_;
    bool hasCharges_;
    bool useSybylTypes_;
    bool prependExt_;
    Mol2File file_; 
};
#endif
