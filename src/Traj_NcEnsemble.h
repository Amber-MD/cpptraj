#ifndef INC_TRAJ_NCENSEMBLE_H
#define INC_TRAJ_NCENSEMBLE_H
#ifdef BINTRAJ
#ifdef ENABLE_SINGLE_ENSEMBLE
#include "TrajectoryIO.h"
#include "NetcdfFile.h"
/// Read/write single NetCDF ensemble trajectory.
class Traj_NcEnsemble : public TrajectoryIO, private NetcdfFile {
  public:
    Traj_NcEnsemble();
    ~Traj_NcEnsemble();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_NcEnsemble(); }
    static void WriteHelp();
    static void ReadHelp();
    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(FileName const&, Topology*);
    int setupTrajout(FileName const&, Topology*, CoordinateInfo const&, int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,Frame&);
    int readVelocity(int, Frame&);
    int readForce(int, Frame&) { return 1; } //FIXME implement
    int writeFrame(int,Frame const&);
    void Info();
    int processWriteArgs(ArgList&, DataSetList const&);
    int processReadArgs(ArgList&);
    bool CanProcessEnsemble() { return true; }
    int readArray(int, FrameArray&);
    int writeArray(int, FramePtrArray const&);
  private:
    float *Coord_;
    FileName filename_;
    int ensembleStart_;
    int ensembleEnd_;
    bool readAccess_;
    bool useVelAsCoords_;
    bool useFrcAsCoords_;
};
#endif
#endif
#endif
