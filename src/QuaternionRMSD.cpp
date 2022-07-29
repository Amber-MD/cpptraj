#include "QuaternionRMSD.h"
#include "qcprot.h"
#include "Constants.h"
#include "CpptrajStdio.h"

/** Calculate the coordinate covariance matrix and perform the
  * quaternion RMSD calculation.
  * \param Ref Reference coordinates centered at the origin; will not be modified.
  * \param Tgt Target coordinates; will be centered at the origin.
  * \param U If specified, will hold rotation vectors in columns.
  * \param Trans Will be set to translation from Target to origin.
  * \param useMass If true, weight by Target atom masses.
  */
double do_quaternion_rmsd(Frame const& Ref, Frame& Tgt,
                          double* U, double* Trans,
                          bool useMass)
{
  // Center Tgt on the origin
  Trans[0] = 0;
  Trans[1] = 0;
  Trans[2] = 0;
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
  Trans[0] = -Trans[0];
  Trans[1] = -Trans[1];
  Trans[2] = -Trans[2];
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

  double rmsd;
  //int err =
  Cpptraj::QCPRot::FastCalcRMSDAndRotation(U, rot.Dptr(), &rmsd, mwss, total_mass);
  //mprintf("DEBUG: qcprot returned %i\n", err);
  return rmsd;
}

/** Calculate quaternion RMSD with translation vector and rotation matrix. */
double QuaternionRMSD_CenteredRef(Frame const& Ref, Frame& Tgt,
                                  Matrix_3x3& U, Vec3& Trans,
                                  bool useMass)
{
  return do_quaternion_rmsd(Ref, Tgt, U.Dptr(), Trans.Dptr(), useMass);
}

/** Calculate quaternion RMSD only. */
double QuaternionRMSD_CenteredRef(Frame const& Ref, Frame& Tgt, bool useMass)
{
  Vec3 Trans;
  return do_quaternion_rmsd(Ref, Tgt, 0, Trans.Dptr(), useMass);
}
