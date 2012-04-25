#ifndef INC_INTEGRATE_H
#define INC_INTEGRATE_H
#include <vector>
// Class: Interpolate
class Interpolate {
  public:
    void cubicSpline_coeff(double *, double *, int);
    int cubicSpline_eval(double*, double*, int);

    void Set_meshX(double, double, int);
    double Integrate_Trapezoid();

    int Mesh_Size() { return mesh_size_; }
    double X(int i) { return mesh_x_[i]; }
    double Y(int i) { return mesh_y_[i]; }
  private:
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> d;

    std::vector<double> mesh_x_;
    std::vector<double> mesh_y_;
    int mesh_size_;
};
#endif
