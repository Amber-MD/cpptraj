#ifndef INC_PARAMETERTYPES_H
#define INC_PARAMETERTYPES_H
class ParmBondType {
  public:
    ParmBondType() : rk_(0), req_(0) {}
    ParmBondType(double rk, double req) : rk_(rk), req_(req) {}
    inline double Rk() { return rk_; }
    inline double Req() { return req_; }
  private:
    double rk_;
    double req_;
};
#endif
