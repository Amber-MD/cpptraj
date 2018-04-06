#ifndef INC_QUATERNIONRMSD_H
#define INC_QUATERNIONRMSD_H
#include "Frame.h"
class QuaternionRMSD {
  public:
    QuaternionRMSD() : Xtgt_(0), Xref_(0), M_(0) {}
    ~QuaternionRMSD();
    void Clear();
    int Init(int, std::vector<double> const&);
    double RMSD_CenteredRef(Frame const&, Frame&, Matrix_3x3&, Vec3&);
  private:
    double** Xtgt_;
    double** Xref_;
    double* M_;
    int len_;
};
#endif
