#ifndef INC_PARAMETERTYPES_H
#define INC_PARAMETERTYPES_H
class ParmBondType {
  public:
    ParmBondType() : rk_(0), req_(0) {}
    ParmBondType(double rk, double req) : rk_(rk), req_(req) {}
    inline const double Rk() const { return rk_; }
    inline const double Req() const { return req_; }
  private:
    double rk_;
    double req_;
};
typedef ParmBondType ParmAngleType;
class ParmNonbondType {
  public:
    ParmNonbondType() : A_(0), B_(0) {}
    ParmNonbondType(double rk, double req) : A_(rk), B_(req) {}
    inline const double A() const { return A_; }
    inline const double B() const { return B_; }
  private:
    double A_;
    double B_;
};
#endif
