#ifndef INC_INTEGRATE_H
#define INC_INTEGRATE_H
// Class: Interpolate
class Interpolate {
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> d;
  public:
    void cubicSpline_coeff(double *, double *, int);
    int cubicSpline_eval(double*,double*,int,double*,double*,int);
};
#endif
