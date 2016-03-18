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

    bool isBigEndian_;   /// True if byte order is reversed
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
    size_t frameSize_;
    size_t headerBytes_;
    size_t arraySize_;
    float* farray_;
    double* darray_;

    void GmxInfo();
    bool IsTRX(CpptrajFile&);
    int read_int(int&);
    int write_int(int);
    int read_real(float&);
    int write_real(float);
    std::string read_string();
    int ReadBox(double*);
    int ReadTrxHeader();
    int ReadAtomVector(double*, int);
    void AllocateCoords();

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
    int readForce(int, Frame&)     { return 1; }  // TODO support this
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
};
#endif
