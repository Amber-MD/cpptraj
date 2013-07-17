#ifndef INC_DATASET_MESH_H
#define INC_DATASET_MESH_H
#include <vector>
#include "DataSet_1D.h"
/// Hold a mesh of X-Y values
class DataSet_Mesh : public DataSet_1D {
  public:
    DataSet_Mesh() : DataSet_1D(XYMESH, 12, 4) {}
    static DataSet* Alloc() { return (DataSet*)new DataSet_Mesh();}
    // ----- DataSet functions -------------------
    size_t Size()            const { return mesh_x_.size();     }
    int Sync()                     { return 0;                  }
    void Info()              const { return;                    }
    // ----- DataSet_1D functions ----------------
    int Allocate1D(size_t);
    void Add( size_t, const void* ) {} // TODO: Implement?
    double Dval(size_t idx)  const { return mesh_y_[idx];       }
    void WriteBuffer(CpptrajFile&, size_t) const;
    // -------------------------------------------
    double X(int i) const { return mesh_x_[i]; }
    double Y(int i) const { return mesh_y_[i]; }
    void CalculateMeshX(int,double,double);
    int SetMeshY(DataSet_1D const&);
    double Integrate_Trapezoid( DataSet_Mesh& );
    double Integrate_Trapezoid();
  private:
    std::vector<double> mesh_x_;
    std::vector<double> mesh_y_;
};
#endif
