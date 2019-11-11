#ifndef INC_TRAJ_TNG_H
#define INC_TRAJ_TNG_H
#ifdef HAS_TNGFILE
#include <tng/tng_io.h>
#include "TrajectoryIO.h"
#include "FileName.h"
/// Read Gromacs TNG trajectories 
class Traj_GmxTng : public TrajectoryIO {
  public:
    Traj_GmxTng();
    ~Traj_GmxTng();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_GmxTng(); }
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

    void convertArray(double*, float*, unsigned int, double) const;
    int getNextBlocks(int64_t&);
    int readValues(int64_t, int64_t&, double&, char&);

    typedef std::vector<int64_t> Iarray;

    tng_trajectory_t traj_;  ///< The TNG trajectory file object
    void* values_;           ///< Temporary array for reading in values from TNG
    int64_t tngatoms_;       ///< Number of atoms in the TNG trajectory file.
    int64_t tngframes_;      ///< Number of *MD sim( frames in the TNG trajectory file.
    int64_t tngsets_;        ///< Number of actual frames in the TNG trajectory file.
    int64_t current_frame_;  ///< The current frame (relative to MD sim, not the trajectory!)
    int current_set_;        ///< The current set in tng file. Determines if we need to seek.
    int64_t next_nblocks_;   ///< The number of data blocks in the next frame
    int64_t* next_blockIDs_; ///< Array containing block IDs in next frame
    double tngfac_;          ///< Coordinates scaling factor
    bool isOpen_;            ///< Calling the TNG library close routine if file is not open is an error, so keep track ourselves.
    FileName filename_;      ///< File name, for openTrajin
    Iarray blockIds_;        ///< Currently active block IDs
};
#endif /* HAS_TNGFILE */
#endif
