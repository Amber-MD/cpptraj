#ifndef INC_TRAJ_TNG_H
#define INC_TRAJ_TNG_H
#ifndef NO_TNGFILE
#include "TrajectoryIO.h"
# include <tng/tng_io.h>
/// Read Gromacs TNG trajectories 
class Traj_TNG : public TrajectoryIO {
  public:
    Traj_TNG();
    ~Traj_TNG();
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

    void convertArray(double*, float*, unsigned int) const;

    tng_trajectory_t traj_; ///< The TNG trajectory file object
    float* ftmp_;           ///< Temporary array for reading in coordinates
    int64_t tngatoms_;      ///< Number of atoms in the TNG trajectory file.
    double tngfac_;         ///< Coordinates scaling factor
    bool isOpen_;           ///< Calling the TNG library close routine if file is not open is an error, so keep track ourselves.
};
#endif /* NO_TNGFILE */
#endif
