#ifndef TRAJIN_SINGLE_H
#define TRAJIN_SINGLE_H
#include "Trajin.h"
/// Class for reading in single trajectories.
class Trajin_Single : public Trajin {
  public:
    Trajin_Single();
    ~Trajin_Single();
    /// Set up trajectory for reading, optionally checking box info.
    int SetupTrajRead(std::string const&, ArgList&, Topology*, bool);
    // ----- Inherited functions -----------------
    /// Set up trajectory for reading, check box info.
    int SetupTrajRead(std::string const&, ArgList&, Topology*);
    int ReadTrajFrame(int, Frame&);
    int BeginTraj(bool);
    void EndTraj();
    void PrintInfo(int) const;
    CoordinateInfo const& TrajCoordInfo() const { return cInfo_; }

    // NOTE: The following are currently for testing Trajin_Ensemble
    void EnsembleInfo() const {} 
    int EnsembleSetup(FrameArray&,FramePtrArray&) {return 1;}
    int ReadEnsemble(int,FrameArray&,FramePtrArray&) {return 1;}
    bool  BadEnsemble() const { return true; }
    // -------------------------------------------
    std::string Title() {
      if (trajio_==0) return std::string("");
      else return trajio_->Title();
    }
  private:
    TrajectoryIO* trajio_; ///< Hold class that will interface with traj format.
    TrajectoryIO* velio_;  ///< Hold class that will interface with opt. mdvel file.
    CoordinateInfo cInfo_; ///< Hold coordinate metadata.
    bool trajIsOpen_;      ///< True if trajectory is open. 
};
#endif
