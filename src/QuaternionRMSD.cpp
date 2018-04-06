#include "QuaternionRMSD.h"
#include "qcprot.h"
#include "Constants.h"
#include "CpptrajStdio.h"

QuaternionRMSD::~QuaternionRMSD() {
  Clear();
}

void QuaternionRMSD::Clear() {
  if (M_ != 0) {
    delete[] M_;
    M_ = 0;
  }
  if (Xtgt_ != 0) {
    for (int i = 0; i < 3; i++) {
      delete[] Xtgt_[i];
      delete[] Xref_[i];
    }
    delete[] Xtgt_;
    delete[] Xref_;
    Xtgt_ = 0;
    Xref_ = 0;
  }
  len_ = 0;
}

int QuaternionRMSD::Init(int natom, std::vector<double> const& mass)
{
  Clear();
  len_ = natom;
  Xtgt_ = new double*[ 3 ];
  Xref_ = new double*[ 3 ];
  for (int i = 0; i < 3; i++) {
    Xtgt_[i] = new double[ len_ ];
    Xref_[i] = new double[ len_ ];
  }
  if (!mass.empty()) {
    M_ = new double[ mass.size() ];
    std::copy(mass.begin(), mass.end(), M_);
  }
  return 0;
}

double QuaternionRMSD::RMSD_CenteredRef(Frame const& Ref, Frame& Tgt,
                                        Matrix_3x3& U, Vec3& Trans)
{
  Trans.Zero();
  double* X_ = Tgt.xAddress();
  const double* R_ = Ref.xAddress();
  int ncoord = len_ * 3;
  double total_mass;
  if (M_ != 0) {
    total_mass = 0.0;
    int im = 0;
    for (int ix = 0; ix < ncoord; ix += 3, im++) {
      double mass = M_[im];
      total_mass += mass;
      Trans[0] += (X_[ix  ] * mass);
      Trans[1] += (X_[ix+1] * mass);
      Trans[2] += (X_[ix+2] * mass);
    }
  } else {
    total_mass = (double)len_;
    int im = 0;
    for (int ix = 0; ix < ncoord; ix += 3, im++) {
      Trans[0] += X_[ix  ];
      Trans[1] += X_[ix+1];
      Trans[2] += X_[ix+2];
     }
  }
  if (total_mass < Constants::SMALL) {
    mprinterr("Error: RMSD: Divide by zero.\n");
    return -1;
  }
  Trans[0] /= total_mass;
  Trans[1] /= total_mass;
  Trans[2] /= total_mass;
  Trans.Neg();
  Tgt.Translate(Trans);
  // Save coordinates
  int ix = 0;
  for (int im = 0; im < len_; im++, ix += 3)
  {
    Xtgt_[0][im] = X_[ix  ];
    Xtgt_[1][im] = X_[ix+1];
    Xtgt_[2][im] = X_[ix+2];
    Xref_[0][im] = R_[ix  ];
    Xref_[1][im] = R_[ix+1];
    Xref_[2][im] = R_[ix+2];
  }

  return CalcRMSDRotationalMatrix( Xref_, Xtgt_, len_, U.Dptr(), M_ );
}
