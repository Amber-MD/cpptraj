#ifndef INC_VEC3_H
#define INC_VEC3_H
/// Designed to hold array of size 3 (like XYZ coord etc).
class Vec3 {
  public:
    Vec3() {
      Zero(); 
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
    Vec3(const double* XYZ) {
      V_[0] = XYZ[0];
      V_[1] = XYZ[1];
      V_[2] = XYZ[2];
    }
    Vec3(const int* XYZ) {
      V_[0] = (double)XYZ[0];
      V_[1] = (double)XYZ[1];
      V_[2] = (double)XYZ[2];
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
    void operator*=(double xIn) {
      V_[0] *= xIn;
      V_[1] *= xIn;
      V_[2] *= xIn;
    }
    void operator+=(double xIn) {
      V_[0] += xIn;
      V_[1] += xIn;
      V_[2] += xIn;
    }
    void operator-=(const Vec3& rhs) {
      V_[0] -= rhs.V_[0];
      V_[1] -= rhs.V_[1];
      V_[2] -= rhs.V_[2];
    }
    void operator+=(const Vec3& rhs) {
      V_[0] += rhs.V_[0];
      V_[1] += rhs.V_[1];
      V_[2] += rhs.V_[2];
    }
    Vec3 operator-(const Vec3& rhs) const {
      Vec3 tmp( *this );
      tmp -= rhs;
      return tmp;
    }
    Vec3 operator+(const Vec3& rhs) const {
      Vec3 tmp( *this );
      tmp += rhs;
      return tmp;
    }
    double operator*(const Vec3& rhs) const { // Dot product
      return ( (V_[0]*rhs.V_[0]) + (V_[1]*rhs.V_[1]) + (V_[2]*rhs.V_[2]) );
    }
    Vec3 operator*(double xIn) {
      Vec3 tmp( *this );
      tmp *= xIn;
      return tmp;
    }
    double operator[](int idx) const { return V_[idx]; }
    double Magnitude2() {
      double x = V_[0] * V_[0];
      double y = V_[1] * V_[1];
      double z = V_[2] * V_[2];
      return (x + y + z);
    }
    void Zero() {
      V_[0] = 0.0;
      V_[1] = 0.0;
      V_[2] = 0.0;
    }
    void Neg() {
      V_[0] = -V_[0];
      V_[1] = -V_[1];
      V_[2] = -V_[2];
    }
    void SetVec(double vx, double vy, double vz) {
      V_[0] = vx;
      V_[1] = vy;
      V_[2] = vz;
    }
    // TODO: Eliminate this routine
    double* Dptr() { return V_; }
  private:
    double V_[3];
};
#endif
