#include "DataSet_Mesh.h"
#include "CpptrajStdio.h"

int DataSet_Mesh::Allocate1D( size_t sizeIn ) {
  mesh_x_.reserve( sizeIn );
  mesh_y_.reserve( sizeIn );
  return 0;
}

void DataSet_Mesh::WriteBuffer(CpptrajFile &cbuffer, size_t frame) const {
  if (frame >= mesh_x_.size())
    cbuffer.Printf(data_format_, 0.0);
  else
    cbuffer.Printf(data_format_, mesh_y_[frame]);
}

void DataSet_Mesh::CalculateMeshX(int sizeIn, double ti, double tf) {
  mesh_x_.resize( sizeIn, 0 );
  mesh_y_.resize( sizeIn, 0 );
  double s = (ti + tf)/2;
  double d = (tf - ti)/2;
  for (int i = 0; i < sizeIn; i++)
    mesh_x_[i] = s + d*((double) (2*i + 1 - sizeIn)/(sizeIn - 1));
}

int DataSet_Mesh::SetMeshY(DataSet_1D const& dsIn) {
  if (dsIn.Size() != mesh_x_.size()) {
    mprintf("Warning: Input data set %s size %zu != mesh %s size %zu\n",
            dsIn.Legend().c_str(), dsIn.Size(), Legend().c_str(), Size());
  }
  size_t maxSize = std::min( dsIn.Size(), Size() );
  for (int i = 0; i < (int)maxSize; i++)
    mesh_y_[i] = dsIn.Dval(i);
  return 0;
}

double DataSet_Mesh::Integrate_Trapezoid( DataSet_Mesh& sumOut ) {
  double sum = 0.0;
  int mesh_size = (int)mesh_x_.size();
  if (mesh_size < 2) return 0.0;
  // Give output data set the same X mesh
  sumOut.mesh_x_ = mesh_x_;
  sumOut.mesh_y_.resize( mesh_x_.size() );
  sumOut.mesh_y_[0] = 0.0;
  for (int i = 1; i < mesh_size; i++) {
      double b_minus_a = (mesh_x_[i] - mesh_x_[i - 1]);
      sum += (b_minus_a * (mesh_y_[i - 1] + mesh_y_[i]) * 0.5);
      sumOut.mesh_y_[i] = sum;
  }
  return sum;
}

double DataSet_Mesh::Integrate_Trapezoid() {
  double sum = 0.0;
  int mesh_size = (int)mesh_x_.size();
  if (mesh_size < 2) return 0.0;
  for (int i = 1; i < mesh_size; i++) {
      double b_minus_a = (mesh_x_[i] - mesh_x_[i - 1]);
      sum += (b_minus_a * (mesh_y_[i - 1] + mesh_y_[i]) * 0.5);
  }
  return sum;
}
