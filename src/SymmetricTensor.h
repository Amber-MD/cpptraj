#ifndef INC_SYMMETRICTENSOR_H
#define INC_SYMMETRICTENSOR_H
/// Class for holding elements of a 3x3 symmetric tensor (6 elements total)
/**  U11 U12 U13
  *  U12 U22 U23
  *  U13 U23 U33
  * stored as
  * U11 U22 U33 U12 U13 U23
  */
template <class T> class SymmetricTensor {
  public:
    /// CONSTRUCTOR
    SymmetricTensor() {}
    /// CONSTRUCTOR - take pointer to array of size 6
    SymmetricTensor(T const* p) {
      U_[0] = p[0];
      U_[1] = p[1];
      U_[2] = p[2];
      U_[3] = p[3];
      U_[4] = p[4];
      U_[5] = p[5];
    }

    T const& U11() const { return U_[0]; }
    T const& U22() const { return U_[1]; }
    T const& U33() const { return U_[2]; }
    T const& U12() const { return U_[3]; }
    T const& U13() const { return U_[4]; }
    T const& U23() const { return U_[5]; }

    T const* Ptr() const { return U_; }
  private:
    T U_[6];
};
#endif
