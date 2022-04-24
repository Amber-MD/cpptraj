#ifndef INC_COORDCOVARMATRIX_H
#define INC_COORDCOVARMATRIX_H
template <class T> class CoordCovarMatrix {
  public:
    CoordCovarMatrix() {}
    /// Calculate mass-weighted covariance matrix of centered Coords and Reference (R = Xt * Ref)
    T CalcCovariance_MassWt(int nselected, const int* imask, T const* Ref, T const* Tgt, T const* Mass) {
      T mwss = 0.0;
      for (unsigned int idx = 0; idx != 9; idx++)
        rot_[idx] = 0;
      for (int aidx = 0; aidx != nselected; aidx++)
      {
        int at = imask[aidx];
        int i = at * 3;
        T xt = Tgt[i  ];
        T yt = Tgt[i+1];
        T zt = Tgt[i+2];
        T xr = Ref[i  ];
        T yr = Ref[i+1];
        T zr = Ref[i+2];
        T atom_mass = Mass[at];
        mwss += atom_mass * ( (xt*xt)+(yt*yt)+(zt*zt)+(xr*xr)+(yr*yr)+(zr*zr) );
        rot_[0] += atom_mass*xt*xr;
        rot_[1] += atom_mass*xt*yr;
        rot_[2] += atom_mass*xt*zr;
        rot_[3] += atom_mass*yt*xr;
        rot_[4] += atom_mass*yt*yr;
        rot_[5] += atom_mass*yt*zr;
        rot_[6] += atom_mass*zt*xr;
        rot_[7] += atom_mass*zt*yr;
        rot_[8] += atom_mass*zt*zr;
      }
      return mwss;
    }
  private:
    T rot_[9]; ///< Hold coordinate covariance matrix, row-major
};
#endif
