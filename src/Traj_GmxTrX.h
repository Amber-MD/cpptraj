#ifndef INC_TRAJ_GMXTRR_H
#define INC_TRAJ_GMXTRR_H
#include "TrajectoryIO.h"
/// Read/write Gromacs TRR/TRJ trajectories
class Traj_GmxTrX : public TrajectoryIO {
  public:
    Traj_GmxTrX();
    ~Traj_GmxTrX();
    static TrajectoryIO* Alloc() { return (TrajectoryIO*)new Traj_GmxTrX(); }
  private:
    enum FormatType { TRR = 0, TRJ };
    static const int Magic_;

    bool isBigEndian_;   /// True if byte order is reversed
    CpptrajFile file_;
    FormatType format_;

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
    float dt_;
    float lambda_;
    size_t frameSize_;
    size_t headerBytes_;
    float* farray_;

    void GmxInfo();
    bool IsTRX(CpptrajFile&);
    int read_int(int&);
    int read_real(float&);
    std::string read_string();
    int ReadBox(double*);
    int ReadTrxHeader();
    int ReadAtomVector(double*, int);

    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(std::string const&, Topology*);
    int setupTrajout(std::string const&, Topology*, int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    void Info();
    int readVelocity(int, double*);
    int processWriteArgs(ArgList&) { return 0; }
    int processReadArgs(ArgList&)  { return 0; }
    int readIndices(int,int*)      { return 1; }
};
#endif
