#ifndef INC_DATASET_MESH_H
#define INC_DATASET_MESH_H
#include <vector>
#include "DataSet_1D.h"
/// Hold a mesh of X-Y values
class DataSet_Mesh : public DataSet_1D {
  public:
    DataSet_Mesh() : DataSet_1D(XYMESH, 12, 4) {}
    /// Construct mesh with preset X values
    DataSet_Mesh(int,double,double);
    static DataSet* Alloc() { return (DataSet*)new DataSet_Mesh();}
    // ----- DataSet functions -------------------
    size_t Size()            const { return mesh_x_.size();     }
    int Sync()                     { return 0;                  }
    void Info()              const { return;                    }
    // ----- DataSet_1D functions ----------------
    int Allocate1D(size_t);
    void Add( size_t, const void* ) {} // TODO: Implement?
    double Dval(size_t idx)  const { return mesh_y_[idx];       }
    double Xcrd(size_t idx)  const { return mesh_x_[idx];       }
    void WriteBuffer(CpptrajFile&, size_t) const;
    // -------------------------------------------
    double X(int i) const { return mesh_x_[i]; }
    double Y(int i) const { return mesh_y_[i]; }
    /// Calculate mesh X values given size, start, and end values.
    void CalculateMeshX(int,double,double);
    /// Set mesh Y values from input data set.
    int SetMeshY(DataSet_1D const&);
    /// Integrate the mesh, compute cumulative sum
    double Integrate_Trapezoid( DataSet_Mesh& ) const;
    /// Integrate the mesh
    double Integrate_Trapezoid() const;
    /// Set mesh with splined values based on input X and Y values.
    int SetSplinedMeshY(std::vector<double> const&, std::vector<double> const&);
    /// Set mesh with splined values based on input DataSet.
    int SetSplinedMesh(DataSet_1D const&);
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
#endif
