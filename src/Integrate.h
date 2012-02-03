#ifndef INC_INTEGRATE_H
#define INC_INTEGRATE_H
//#include <vector>
// Class: Interpolate
/*class Interpolate {
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> d;
  public:*/
    void cubicSpline_coeff(double *, double *, int, double*, double*, double*);
    int cubicSpline_eval(double*,double*,int,double*,double*,double*,double*,double*,int);
//};
    double integrate_trapezoid(double *, double *, int);
#endif
