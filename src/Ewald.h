#ifndef INC_EWALD_H
#define INC_EWALD_H
class Ewald {
  public:
    Ewald();
    void FindEwaldCoefficient(double,double);
  private:
    double erfc_func(double) const;

    double sumq_;
    double sumq2_;
    double ew_coeff_;
};
#endif
