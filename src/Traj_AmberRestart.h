#ifndef INC_TRAJ_AMBERRESTART_H
#define INC_TRAJ_AMBERRESTART_H
#include "TrajectoryIO.h"
#include "BufferedFile.h"
// Class: Traj_AmberRestart.h
/// Reads and writes formatted (ASCII text) amber
class Traj_AmberRestart : public TrajectoryIO {
  public:
    Traj_AmberRestart();
    static TrajectoryIO* Alloc() { return (TrajectoryIO*)new Traj_AmberRestart(); }
    // AmberRestart-specific functions
    void SetNoVelocity();
  private:
    int debug_;
    int restartAtoms_;     ///< Number of atoms in restart file
    int natom3_;           ///< Number of coords
    int numBoxCoords_;     ///< Number of box coords (3 or 6)
    size_t coordSize_;     ///< Size of coords in bytes, for reading past coords.
    double restartTime_;   ///< Time in restart file, read in
    double restartTemp_;   ///< (Optional) replica temperature, read in.
    double time0_;         ///< For writes, restart time offset
    double dt_;            ///< For writes, restart timestep (scaling)
    bool singleWrite_;     ///< If false, frame # will be appended to output filename
    bool HasT_;
    bool HasV_;
    BufferedFile file_;
    std::string title_;

    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(std::string const&, Topology*, TrajInfo&);
    int setupTrajout(std::string const&, Topology*, int, TrajInfo const&,bool);
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int readVelocity(int, double*);
    int writeFrame(int,double*,double*,double*,double);
    int processWriteArgs(ArgList&);
    void info();

    int getBoxAngles(const char*, Box&);
    int readIndices(int,int*) { return 1; }
    int processReadArgs(ArgList&) { return 0; }
};
#endif
