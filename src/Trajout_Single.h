#ifndef INC_TRAJOUT_SINGLE_H
#define INC_TRAJOUT_SINGLE_H
#include "OutputTrajCommon.h"
/// Write out 1 frame at a time to a single file.
/** Note that unlike Trajin, there is really no point in having
  * a single frame written to multiple files since this is handled
  * by TrajoutList. Therefore this class doesnt need to inherit
  * from a base class currently. The Trajout_Single name is used
  * however in case we do want to inherit in the future in a manner
  * similar to Trajin.
  */
class Trajout_Single {
  public:
    Trajout_Single() : trajio_(0), debug_(0) {}
    ~Trajout_Single();
    void SetDebug(int d) { debug_ = d; }
    // ----- Inherited functions -----------------
    /// Prepare trajectory for writing to the given format, but no Topology setup.
    int InitTrajWrite(std::string const&, ArgList const&, TrajectoryFile::TrajFormatType);
    /// Peform Topology-related setup for trajectory and open. TODO const&
    int SetupTrajWrite(Topology*, CoordinateInfo const&, int);
    /// Close output trajectory.
    void EndTraj();
    /// Write a single frame.
    int WriteSingle(int, Frame const&);
    /// Print information on trajectory to be written.
    void PrintInfo(int) const;
    // -------------------------------------------
    OutputTrajCommon Traj() const { return traj_; }
    /// Init and setup/open traj.
    int PrepareTrajWrite(std::string const&, ArgList const&, Topology*,
                         CoordinateInfo const&, int, TrajectoryFile::TrajFormatType);
    /// Init and setup/open traj for writing to STDOUT (e.g. ambpdb mode)
    int PrepareStdoutTrajWrite(ArgList const&, Topology*, CoordinateInfo const&, int,
                               TrajectoryFile::TrajFormatType);
    /// Init traj; if given, append ensemble number to name
    int InitEnsembleTrajWrite(std::string const&, ArgList const&,
                              TrajectoryFile::TrajFormatType, int);
    /// Init and setup/open traj; if given, append ensemble number to name
    int PrepareEnsembleTrajWrite(std::string const&, ArgList const&, Topology*,
                                 CoordinateInfo const&, int,
                                 TrajectoryFile::TrajFormatType, int);
  private:
    int InitTrajout(std::string const&, ArgList const&, TrajectoryFile::TrajFormatType);

    OutputTrajCommon traj_;
    TrajectoryIO* trajio_;
    int debug_;
};
#endif
