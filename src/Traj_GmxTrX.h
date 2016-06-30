#ifndef INC_TRAJ_GMXTRR_H
#define INC_TRAJ_GMXTRR_H
#include "TrajectoryIO.h"
/// Read/write Gromacs TRR/TRJ trajectories
class Traj_GmxTrX : public TrajectoryIO {
  public:
    Traj_GmxTrX();
    ~Traj_GmxTrX();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_GmxTrX(); }
    static void WriteHelp();
  private:
    enum FormatType { TRR = 0, TRJ };
    static const int Magic_;

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
    int processWriteArgs(ArgList&);
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

    void GmxInfo();
    int DetermineEndian(int);
    bool IsTRX(CpptrajFile&);
    int read_int(int&);
    int write_int(int);
    int read_real(float&);
    int write_real(float);
    std::string read_string();
    int ReadBox(double*);
    int ReadTrxHeader(int&);
    void AllocateCoords();

    bool swapBytes_;   ///< True if byte order needs to be reversed
    bool isBigEndian_; ///< True if file is big-endian.
    CpptrajFile file_;
    FormatType format_;

    double dt_;
    int ir_size_;
    int e_size_;
    int box_size_;
    int vir_size_;
    int pres_size_;
    int top_size_;
    int sym_size_;
    int x_size_;
    int v_size_;
    int f_size_;
    int natoms_;
    int natom3_;
    int step_;
    int nre_;
    int precision_;
    float timestep_;
    float lambda_;
    size_t frameSize_;   ///< Size of single trajectory frame in bytes
    size_t headerBytes_; ///< Size of header in bytes
    size_t arraySize_;   ///< # elements in {d|f}array_; total # of position/veloc/force coords.
    float* farray_;      ///< Array for reading/writing single precision.
    double* darray_;     ///< Array for reading/writine double precision.

    static const double GMX_FRC_TO_AMBER;
    static const double AMBER_FRC_TO_GMX;
    static const double GMX_VEL_TO_AMBER;
    static const double AMBER_VEL_TO_GMX;
};
#endif
