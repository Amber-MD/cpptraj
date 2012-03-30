#ifndef INC_VEC3_H
#define INC_VEC3_H
/// Vector of length 3
class Vec3 {
  public:
    Vec3() : x_(0), y_(0), z_(0) {}
    Vec3(double *XYZ) : x_(XYZ[0]), y_(XYZ[1]), z_(XYZ[2]) {}
    Vec3(const Vec3& rhs) : x_(rhs.x_), y_(rhs.y_), z_(rhs.z_) {}
    Vec3 &operator=(const Vec3& rhs) {
      if (&rhs==this) return *this;
      x_ = rhs.x_;
      y_ = rhs.y_;
      z_ = rhs.z_;
      return *this;
    }
    Vec3 &operator-=(const Vec3& rhs) {
      x_ -= rhs.x_;
      y_ -= rhs.y_;
      z_ -= rhs.z_;
      return *this;
    }
    double LengthSquared() {
      return ( (x_*x_)+(y_*y_)+(z_*z_) );
    }
  private:
    double x_;
    double y_;
    double z_;
};
#endif
