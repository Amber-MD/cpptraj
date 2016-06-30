#ifndef INC_TRAJ_GMXXTC_H
#define INC_TRAJ_GMXXTC_H
#include "TrajectoryIO.h"
#ifndef NO_XDRFILE
# include <xdrfile_xtc.h>
#endif
/// Read/write Gromacs XTC trajectories
class Traj_GmxXtc : public TrajectoryIO {
  public:
    Traj_GmxXtc();
    ~Traj_GmxXtc();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_GmxXtc(); }
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
    int readVelocity(int, Frame&);
    int readForce(int, Frame&);
    int processWriteArgs(ArgList&) { return 0; }
    int processReadArgs(ArgList&)  { return 0; }
#   ifdef MPI
    // Parallel functions
    int parallelOpenTrajin(Parallel::Comm const&);
    int parallelOpenTrajout(Parallel::Comm const&);
    int parallelSetupTrajout(FileName const&, Topology*, CoordinateInfo const&,
                             int, bool, Parallel::Comm const&);
    int parallelReadFrame(int, Frame&);
    int parallelWriteFrame(int, Frame const&);
    void parallelCloseTraj();
#   endif
#   ifndef NO_XDRFILE
    XDRFILE* xd_; ///< Hold XDR file metadata
    rvec* vec_;   ///< Temporary location for holding XDR frame data
    matrix box_;  ///< Temporary location for holding XDR box data
    int natoms_;  ///< Number of atoms in xdr file
#   endif
};
#endif
