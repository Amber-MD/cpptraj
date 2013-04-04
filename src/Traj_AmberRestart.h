#ifndef INC_TRAJ_AMBERRESTART_H
#define INC_TRAJ_AMBERRESTART_H
#include "TrajectoryIO.h"
#include "BufferedFrame.h"
// Class: Traj_AmberRestart.h
/// Reads and writes formatted (ASCII text) amber
class Traj_AmberRestart : public TrajectoryIO {
  public:
    Traj_AmberRestart();
    static TrajectoryIO* Alloc() { return (TrajectoryIO*)new Traj_AmberRestart(); }
  private:
    int restartAtoms_;     ///< Number of atoms in restart file
    int natom3_;           ///< Number of coords
    int numBoxCoords_;     ///< Number of box coords (3 or 6)
    size_t coordSize_;     ///< Size of coords in bytes, for reading past coords.
    double restartTime_;   ///< Time in restart file, read in
    double restartTemp_;   ///< (Optional) replica temperature, read in.
    double time0_;         ///< For writes, restart time offset
    double dt_;            ///< For writes, restart timestep (scaling)
    bool singleWrite_;     ///< If false, frame # will be appended to output filename
    bool readAccess_;      ///< If true, presence/absence of velocity info is known
    BufferedFrame file_;

    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(std::string const&, Topology*);
    int setupTrajout(std::string const&, Topology*, int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int readVelocity(int, double*);
    int writeFrame(int,double*,double*,double*,double);
    int processWriteArgs(ArgList&);
    void Info();

    int getBoxAngles(std::string const&, Box&);
    int readIndices(int,int*) { return 1; }
    int processReadArgs(ArgList&) { return 0; }
};
#endif
