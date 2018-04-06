#include "QuaternionRMSD.h"
#include "qcprot.h"
#include "Constants.h"
#include "CpptrajStdio.h"

double QuaternionRMSD_CenteredRef(Frame const& Ref, Frame& Tgt,
                                  Matrix_3x3& U, Vec3& Trans,
                                  bool useMass, double minScore)
{
  Trans.Zero();
  int ncoord = Ref.size();
  double total_mass;
  if (useMass) {
    total_mass = 0.0;
    int im = 0;
    for (int ix = 0; ix < ncoord; ix += 3, im++) {
      double mass = Tgt.Mass(im);
      total_mass += mass;
      Trans[0] += (Tgt[ix  ] * mass);
      Trans[1] += (Tgt[ix+1] * mass);
      Trans[2] += (Tgt[ix+2] * mass);
    }
  } else {
    total_mass = (double)Ref.Natom();
    int im = 0;
    for (int ix = 0; ix < ncoord; ix += 3, im++) {
      Trans[0] += Tgt[ix  ];
      Trans[1] += Tgt[ix+1];
      Trans[2] += Tgt[ix+2];
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

  // Calculate covariance matrix of Coords and Reference (R = Xt * Ref)
  // Calculate the Kabsch matrix: R = (rij) = Sum(yni*xnj)
  double mwss = 0.0;
  Matrix_3x3 rot(0.0);
  if (useMass) {
    int im = 0;
    for (int i = 0; i < ncoord; i += 3, im++)
    {
      double xt = Tgt[i  ];
      double yt = Tgt[i+1];
      double zt = Tgt[i+2];
      double xr = Ref[i  ];
      double yr = Ref[i+1];
      double zr = Ref[i+2];
      double atom_mass = Tgt.Mass(im);
      mwss += atom_mass * ( (xt*xt)+(yt*yt)+(zt*zt)+(xr*xr)+(yr*yr)+(zr*zr) );
      rot[0] += atom_mass*xt*xr;
      rot[1] += atom_mass*xt*yr;
      rot[2] += atom_mass*xt*zr;
      rot[3] += atom_mass*yt*xr;
      rot[4] += atom_mass*yt*yr;
      rot[5] += atom_mass*yt*zr;
      rot[6] += atom_mass*zt*xr;
      rot[7] += atom_mass*zt*yr;
      rot[8] += atom_mass*zt*zr;
    }
  } else {
    for (int i = 0; i < ncoord; i += 3)
    {
      double xt = Tgt[i  ];
      double yt = Tgt[i+1];
      double zt = Tgt[i+2];
      double xr = Ref[i  ];
      double yr = Ref[i+1];
      double zr = Ref[i+2];
      mwss += ( (xt*xt)+(yt*yt)+(zt*zt)+(xr*xr)+(yr*yr)+(zr*zr) );
      rot[0] += xt*xr;
      rot[1] += xt*yr;
      rot[2] += xt*zr;
      rot[3] += yt*xr;
      rot[4] += yt*yr;
      rot[5] += yt*zr;
      rot[6] += zt*xr;
      rot[7] += zt*yr;
      rot[8] += zt*zr;
    }
  }
  mwss *= 0.5;    // E0 = 0.5*Sum(xn^2+yn^2) 
/*
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
*/
  double rmsd;
  //int err =
  FastCalcRMSDAndRotation(U.Dptr(), rot.Dptr(), &rmsd, mwss, total_mass, minScore);
  //mprintf("DEBUG: qcprot returned %i\n", err);
  return rmsd;
}

double QuaternionRMSD_CenteredRef(Frame const& Ref, Frame& Tgt,
                                  Matrix_3x3& U, Vec3& Trans,
                                  bool useMass)
{
  return QuaternionRMSD_CenteredRef(Ref, Tgt, U, Trans, useMass, -1);
}

double QuaternionRMSD_CenteredRef(Frame const& Ref, Frame& Tgt, bool useMass)
{
  Matrix_3x3 U;
  Vec3 Trans;
  return QuaternionRMSD_CenteredRef(Ref, Tgt, U, Trans, useMass, 99999);
}
