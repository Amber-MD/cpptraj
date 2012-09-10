#ifndef INC_VEC3_H
#define INC_VEC3_H
class Vec3 {
  public:
    Vec3() { 
      V_[0] = 0;
      V_[1] = 0;
      V_[2] = 0;
    }
    Vec3(const Vec3& rhs) {
      V_[0] = rhs.V_[0];
      V_[1] = rhs.V_[1];
      V_[2] = rhs.V_[2];
    }
    Vec3(double vx, double vy, double vz) {
      V_[0] = vx;
      V_[1] = vy;
      V_[2] = vz;
    }
    Vec3& operator=(const Vec3& rhs) {
      if (this == &rhs) return *this;
      V_[0] = rhs.V_[0];
      V_[1] = rhs.V_[1];
      V_[2] = rhs.V_[2];
      return *this;
    }
    void operator/=(double xIn) {
      V_[0] /= xIn;
      V_[1] /= xIn;
      V_[2] /= xIn;
    }
    void operator-=(const Vec3& rhs) {
      V_[0] -= rhs.V_[0];
      V_[1] -= rhs.V_[1];
      V_[2] -= rhs.V_[2];
    }
    double Magnitude() {
      double x = V_[0] * V_[0];
      double y = V_[1] * V_[1];
      double z = V_[2] * V_[2];
      return (x + y + z);
    }
    // TODO: Eliminate this routine
    double* Dptr() { return V_; }
  private:
    double V_[3];
};
#endif
