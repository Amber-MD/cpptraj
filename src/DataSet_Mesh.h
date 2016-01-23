#ifndef INC_DATASET_MESH_H
#define INC_DATASET_MESH_H
#include <vector>
#include "DataSet_1D.h"
/// Hold a mesh of X-Y values
class DataSet_Mesh : public DataSet_1D {
  public:
    DataSet_Mesh() : DataSet_1D(XYMESH, TextFormat(TextFormat::DOUBLE, 12, 4)) {}
    /// Construct mesh with preset X values
    DataSet_Mesh(int,double,double);
    static DataSet* Alloc() { return (DataSet*)new DataSet_Mesh();            }
    void Resize(size_t n)   { mesh_x_.resize(n, 0.0); mesh_y_.resize(n, 0.0); }
    // ----- DataSet functions -------------------
    size_t Size()                       const { return mesh_x_.size();     }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&);
#   endif
    void Info()                         const { return;                    }
    int Allocate(SizeArray const&);
    void Add( size_t, const void* );
    void WriteBuffer(CpptrajFile&, SizeArray const&) const;
    int Append(DataSet*);
    double Coord(unsigned int d, size_t p) const { return mesh_x_[p]; }
    // ----- DataSet_1D functions ----------------
    double Dval(size_t idx)  const { return mesh_y_[idx];       }
    double Xcrd(size_t idx)  const { return mesh_x_[idx];       }
    // -------------------------------------------
    inline void AddXY(double,double);
    double X(int i) const { return mesh_x_[i]; }
    double Y(int i) const { return mesh_y_[i]; }
    /// Set mesh Y value at given index.
    void SetY(int i, double y) { mesh_y_[i] = y; }
    /// Calculate mesh X values given size, start, and end values.
    void CalculateMeshX(int,double,double);
    /// Set mesh X and Y values from input data set.
    int SetMeshXY(DataSet_1D const&); // TODO remove
    /// Set mesh X and Y values from input arrays.
    inline int SetMeshXY(std::vector<double> const&, std::vector<double> const&); // TODO remove
    // -------------------------------------------
    /// Integrate the mesh, compute cumulative sum
    double Integrate_Trapezoid( DataSet_Mesh& ) const;
    /// Integrate the mesh
    double Integrate_Trapezoid() const;
    // -------------------------------------------
    /// Set mesh with splined values based on input X and Y values.
    int SetSplinedMeshY(std::vector<double> const&, std::vector<double> const&);
    /// Set mesh with splined values based on input DataSet.
    int SetSplinedMesh(DataSet_1D const&);
    // -------------------------------------------
    /// Calculate single exponential regression via log and linear regression.
    int SingleExpRegression(double&, double&, double&, CpptrajFile*);
  private:
    void cubicSpline_coeff(std::vector<double> const&, std::vector<double> const&);
    void cubicSpline_eval(std::vector<double> const&, std::vector<double> const&);

    std::vector<double> mesh_x_;
    std::vector<double> mesh_y_;
    // Cubic spline coefficients.
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> d;
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
void DataSet_Mesh::AddXY(double x, double y) {
  mesh_x_.push_back( x );
  mesh_y_.push_back( y );
}
int DataSet_Mesh::SetMeshXY(std::vector<double> const& X, std::vector<double> const& Y) {
  if (mesh_x_.size() != mesh_y_.size()) return 1;
  mesh_x_ = X;
  mesh_y_ = Y;
  return 0;
}
#endif
