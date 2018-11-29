#ifndef INC_TRAJ_XYZ_H
#define INC_TRAJ_XYZ_H
#include "TrajectoryIO.h"
#include "BufferedLine.h"
/// Read simple XYZ trajectories. 
class Traj_XYZ : public TrajectoryIO {
  public:
    Traj_XYZ();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_XYZ(); }
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
    int processWriteArgs(ArgList&);
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

    enum Type { UNKNOWN=0, XYZ, ATOM_XYZ };
    enum TitleType { NO_TITLE = 0, SINGLE, MULTIPLE };

    Type DetermineFormat(std::string&, std::string const&) const;
    inline void ReadTitle();

    static const char* FMT_XYZ_;
    static const char* FMT_ATOM_XYZ_;

    BufferedLine infile_;
    TitleType titleType_;
    int set_;
    const char* fmt_; ///< Format for reading
};
#endif
