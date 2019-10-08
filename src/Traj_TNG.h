#ifndef INC_TRAJ_TNG_H
#define INC_TRAJ_TNG_H
#ifndef NO_TNGFILE
#include "TrajectoryIO.h"
# include <tng/tng_io.h>
/// Read Gromacs TNG trajectories 
class Traj_TNG : public TrajectoryIO {
  public:
    Traj_TNG();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_TNG(); }
    static void WriteHelp();
    static void ReadHelp();
  private:
    // ----- Inherited functions -----------------
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(FileName const&, Topology*);
    int setupTrajout(FileName const&, Topology*, CoordinateInfo const&,int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,Frame&);
    int writeFrame(int,Frame const&);
    void Info();
    int readVelocity(int, Frame&);
    int readForce(int, Frame&);
    int processWriteArgs(ArgList&, DataSetList const&);
    int processReadArgs(ArgList&);
    // -------------------------------------------
#   ifdef MPI
    // ----- Parallel functions ------------------
    int parallelOpenTrajin(Parallel::Comm const&);
    int parallelOpenTrajout(Parallel::Comm const&);
    int parallelSetupTrajout(FileName const&, Topology*, CoordinateInfo const&,
                             int, bool, Parallel::Comm const&);
    int parallelReadFrame(int, Frame&);
    int parallelWriteFrame(int, Frame const&);
    void parallelCloseTraj();
    // -------------------------------------------
#   endif

    tng_trajectory_t traj_; ///< The TNG trajectory file object
    int64_t tngatoms_;      ///< Number of atoms in the TNG trajectory file.
    double tngfac_;         ///< Coordinates scaling factor
};
#endif /* NO_TNGFILE */
#endif
