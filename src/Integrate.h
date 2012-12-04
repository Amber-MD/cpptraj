#ifndef INC_INTEGRATE_H
#define INC_INTEGRATE_H
#include <vector>
// Class: Interpolate
class Interpolate {
  public:
    Interpolate() : mesh_size_(0) {}
    Interpolate(double,double,int);
    void SetMesh_X(double, double, int);
    int SetMesh_Y(std::vector<double> const&, std::vector<double> const&);
    double Integrate_Trapezoid();

    int Mesh_Size() { return mesh_size_; }
    double X(int i) { return mesh_x_[i]; }
    double Y(int i) { return mesh_y_[i]; }
  private:
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> d;
    int mesh_size_;
    std::vector<double> mesh_x_;
    std::vector<double> mesh_y_;

    void CalcMeshX(double,double);
    void cubicSpline_coeff(std::vector<double> const&, std::vector<double> const&);
    void cubicSpline_eval(std::vector<double> const&, std::vector<double> const&);
};
#endif
