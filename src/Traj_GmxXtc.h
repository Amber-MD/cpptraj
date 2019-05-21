#ifndef INC_TRAJ_GMXXTC_H
#define INC_TRAJ_GMXXTC_H
#include "TrajectoryIO.h"
#ifndef NO_XDRFILE
# include <xdrfile_xtc.h>
# include "FileName.h"
#endif
/// Read/write Gromacs XTC trajectories
class Traj_GmxXtc : public TrajectoryIO {
  public:
    Traj_GmxXtc();
    ~Traj_GmxXtc();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_GmxXtc(); }
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
    int readVelocity(int, Frame&);
    int readForce(int, Frame&);
    int processWriteArgs(ArgList&, DataSetList const&);
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
    std::vector<off_t> frameOffsets_; ///< Frame offsets for reading
    XDRFILE* xd_; ///< Hold XDR file metadata
    rvec* vec_;   ///< Temporary location for holding XDR frame data
    matrix box_;  ///< Temporary location for holding XDR box data
    double dt_;   ///< Time step between frames for output
    int natoms_;  ///< Number of atoms in xdr file
    FileName fname_; ///< File name
    float prec_;  ///< Precision
#   endif
};
#endif
