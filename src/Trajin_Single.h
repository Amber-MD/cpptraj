#ifndef INC_TRAJIN_SINGLE_H
#define INC_TRAJIN_SINGLE_H
#include "Trajin.h"
#include "TrajectoryIO.h"
/// Read in 1 frame at a time from a single file.
class Trajin_Single : public Trajin {
  public:
    Trajin_Single();
    ~Trajin_Single();
    // ----- Inherited functions ------------------
    /// Set up trajectory for reading.
    int SetupTrajRead(FileName const&, ArgList&, Topology*);
    /// Read specified frame #.
    int ReadTrajFrame(int, Frame&);
    /// Prepare trajectory for reading.
    int BeginTraj();
    /// Close trajectory.
    void EndTraj();
    /// Print trajectory information.
    void PrintInfo(int) const;
    /// \return trajectory metadata.
    CoordinateInfo const& TrajCoordInfo() const { return cInfo_; }
    // ---------------------------------------------
    std::string const& Title() const { return trajio_->Title(); } //TODO Check for segfault
#   ifdef MPI
    int ParallelBeginTraj( Parallel::Comm const& );
    int ParallelReadTrajFrame(int, Frame&);
    void ParallelEndTraj();
#   endif
  private:
    TrajectoryIO* SetupSeparateTraj(FileName const&, const char*) const;

    TrajectoryIO* trajio_; ///< Hold class that will interface with traj format.
    TrajectoryIO* velio_;  ///< Hold class that will interface with opt. mdvel file.
    TrajectoryIO* frcio_;  ///< Hold class taht will interface with opt. mdfrc file.
    CoordinateInfo cInfo_; ///< Hold coordinate metadata.
};
#endif
