#ifndef INC_CPPTRAJ_DCALC_H
#define INC_CPPTRAJ_DCALC_H
#include "Frame.h" // will bring in Vec3, Matrix_3x3, Box
#include <limits> // std::numeric_limits
#include <cmath> // std::fabs
namespace Cpptraj {
/** Class for implementing multiple distance calcs for a single frame. */
class Dcalc {
  /// Bring each fractional coordinate back into the main unit cell (0 to 1)
  static inline void wrap_frac(Vec3& frac) {
    for (int idx = 0; idx != 3; idx++) {
      while (frac[idx] < 0.0)
        frac[idx] += 1.0;
      while (frac[idx] > 1.0)
        frac[idx] -= 1.0;
    }
  }
  public:
    /// CONSTRUCTOR
    Dcalc( Frame const& frmIn, ImageOption::Type itype) :
      ucell_( frmIn.BoxCrd().UnitCell() ),
      fcell_( frmIn.BoxCrd().FracCell() ),
      boxX_( frmIn.BoxCrd().Param(Box::X) ),
      boxY_( frmIn.BoxCrd().Param(Box::Y) ),
      boxZ_( frmIn.BoxCrd().Param(Box::Z) ),
      xyz1_(0),
      itype_(ImageOption::NO_IMAGE)
    {
      if (boxX_ != 0.0 && boxY_ != 0.0 && boxZ_ != 0.0) {
        itype_ = itype;
      }
    }
    /// Set xyz1 - typically the coordinate from an outer loop
    void SetXYZ1(const double* xyz1) {
      xyz1_ = xyz1;
      if (itype_ == ImageOption::NONORTHO) {
        f1_ = fcell_ * Vec3(xyz1_);
        wrap_frac(f1_);
      }
    }
  private:
    /// \return Dist^2 for xyz2 from previously set xyz1
    double Dist2_ortho(const double* xyz2) const
    {
      // -----------------------------------------
//      if (isOrtho_) {
        double x = xyz1_[0] - xyz2[0];
        double y = xyz1_[1] - xyz2[1];
        double z = xyz1_[2] - xyz2[2];
        x = std::fabs(x);
        y = std::fabs(y);
        z = std::fabs(z);
        while (x > boxX_) x -= boxX_;
        while (y > boxY_) y -= boxY_;
        while (z > boxZ_) z -= boxZ_;
        double D = boxX_ - x;
        if (D < x) x = D;
        D = boxY_ - y;
        if (D < y) y = D;
        D = boxZ_ - z;
        if (D < z) z = D;
        return (x*x + y*y + z*z);
      // -----------------------------------------
    }
    double Dist2_nonortho(const double* xyz2) const
    {
//      } else if (isNonortho_) {
        //return DIST2_ImageNonOrtho(xyz1, xyz2, box.UnitCell(), box.FracCell());
        // Wrap each coordinate in fractional coordinates back to the main unit cell
        // (so each frac coord is between 0 and 1).
        Vec3 f2 = fcell_ * Vec3(xyz2);
        wrap_frac(f2);
        // Calculate f2 to f1 in fractional space
        Vec3 fdelta = f2 - f1_;
        //fdelta.Print("fdelta");
        // Set the min and max cell indices we need to look for based
        // on the f1->f2 vector.
        int imin[3];
        int imax[3];
        for (int idx = 0; idx != 3; idx++) {
          if (fdelta[idx] < 0) {
            imin[idx] = -1;
            imax[idx] = 1;
          } else {
            imin[idx] = 0;
            imax[idx] = 2;
          }
        }
        // Calculate wrapped f2 back in Cartesian space
        Vec3 c2 = ucell_.TransposeMult( f2 );
        // Searching only images we are close to.
        // Check each image of f1 to determine which distance is the smallest.
        double minD2 = std::numeric_limits<double>::max();
        //ixyz[0] = 0;
        //ixyz[1] = 0;
        //ixyz[2] = 0;
        for (int ix = imin[0]; ix != imax[0]; ix++) {
          for (int iy = imin[1]; iy != imax[1]; iy++) {
            for (int iz = imin[2]; iz != imax[2]; iz++) {
              // Image in frac space
              Vec3 f1img = f1_ + Vec3(ix, iy, iz);
              // Image in Cartesian space
              Vec3 c1img = ucell_.TransposeMult( f1img );
              // Calculate distance in Cartesian space
              double dx = (c1img[0] - c2[0]);
              double dy = (c1img[1] - c2[1]);
              double dz = (c1img[2] - c2[2]);
              double D2 = (dx*dx) + (dy*dy) + (dz*dz);
              if (D2 < minD2) {
                minD2 = D2;
                //ixyz[0] = ix;
                //ixyz[1] = iy;
                //ixyz[2] = iz;
              }
            } // END loop over Z image indices
          } // END loop over Y image indices
        } // END loop over X image indices
        //mprintf("minD %12.4f\n", sqrt(minD2));
        return minD2;
    }
    double Dist2_noimage(const double* xyz2) const
    {
//      } else {
      // -----------------------------------------
        //return DIST2_NoImage(xyz1_, xyz2);
        double x = xyz1_[0] - xyz2[0];
        double y = xyz1_[1] - xyz2[1];
        double z = xyz1_[2] - xyz2[2];
        //double D = x*x + y*y + z*z;
        //fprintf(stdout,"Mask1=%8.3f %8.3f %8.3f Mask2=%8.3f %8.3f %8.3f D=%8.3f\n",
        //        a1[0],a1[1],a1[2],a2[0],a2[1],a2[2],D);
        return (x*x + y*y + z*z);
     }
  public:
    double Dist2(const double* xyz2) const {
      switch (itype_) {
        case ImageOption::ORTHO : return Dist2_ortho(xyz2);
        case ImageOption::NONORTHO : return Dist2_nonortho(xyz2);
        case ImageOption::NO_IMAGE : return Dist2_noimage(xyz2);
      }
      return -1.0;
    }
  private:
    Matrix_3x3 const& ucell_;
    Matrix_3x3 const& fcell_;
    Vec3 f1_;
    double boxX_, boxY_, boxZ_;
    const double* xyz1_;
    ImageOption::Type itype_;
};
}
#endif
