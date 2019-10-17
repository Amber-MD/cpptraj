#ifndef INC_TRAJ_GMXDUMP_H
#define INC_TRAJ_GMXDUMP_H
#include "TrajectoryIO.h"
#include "CpptrajFile.h"
/// Can be used to compare output to 'gmx dump'. Currently only for debug.
class Traj_GmxDump : public TrajectoryIO {
  public:
    Traj_GmxDump();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_GmxDump(); }
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

    void indent(int);
    void writeVectorArray(const double*, const char*, int, int, int, int, double);

    CpptrajFile file_;
    int natoms_; ///< Number of atoms in file
    const char* outfmt_; ///< Hold output write format
    bool longFormat_;    ///< If true use the longer format
    bool tngfmt_;        ///< If true use TNG style output format
};
#endif
