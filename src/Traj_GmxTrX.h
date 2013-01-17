#ifndef INC_TRAJ_GMXTRR_H
#define INC_TRAJ_GMXTRR_H
#include "TrajectoryIO.h"
/// Read/write Gromacs TRR/TRJ trajectories
class Traj_GmxTrX : public TrajectoryIO {
  public:
    Traj_GmxTrX();
    static TrajectoryIO* Alloc() { return (TrajectoryIO*)new Traj_GmxTrX(); }
  private:
    enum FormatType { TRR = 0, TRJ };
    static const int Magic_;

    bool isBigEndian_;   /// True if byte order is reversed
    CpptrajFile file_;
    FormatType format_;
    static const int BUF_SIZE = 128;
    char linebuffer_[BUF_SIZE];

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
    int step_;
    int nre_;
    int precision_;
    float dt_;
    float lambda_;

    void GmxInfo();
    bool IsTRX(CpptrajFile&);
    int read_int(int&);
    int read_real(float&);
    std::string read_string();

    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(std::string const&, Topology*);
    int setupTrajout(std::string const&, Topology*, int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    void Info();
    int readVelocity(int, double*) { return 1; }
    int processWriteArgs(ArgList&) { return 0; }
    int processReadArgs(ArgList&)  { return 0; }
    int readIndices(int,int*)      { return 1; }
};
#endif
