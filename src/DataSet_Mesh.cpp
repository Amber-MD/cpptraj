#include <cmath> // regression
#include "DataSet_Mesh.h"
#include "CpptrajStdio.h"
#include "Constants.h" // regression

/// For VoidPtr
static double buf_[2];

/// CONSTRUCTOR - Create X mesh
DataSet_Mesh::DataSet_Mesh(int sizeIn, double ti, double tf) :
  DataSet_1D(XYMESH, TextFormat(TextFormat::DOUBLE, 12, 4))
{
  CalculateMeshX(sizeIn, ti, tf);
}

const void* DataSet_Mesh::VoidPtr(size_t idx) const {
  buf_[0] = mesh_x_[idx];
  buf_[1] = mesh_y_[idx];
  return (const void*)buf_;
}

// DataSet_Mesh::Allocate()
int DataSet_Mesh::Allocate( SizeArray const& sizeIn ) {
  if (!sizeIn.empty()) {
    mesh_x_.reserve( sizeIn[0] );
    mesh_y_.reserve( sizeIn[0] );
  }
  return 0;
}

/** Insert data vIn at frame. */
void DataSet_Mesh::Add(size_t frame, const void* vIn) {
  if (frame > mesh_x_.size()) {
    mesh_x_.resize( frame, 0.0 );
    mesh_y_.resize( frame, 0.0 );
  }
  // Always insert at the end
  // NOTE: No check for duplicate frame values.
  const double* ptr = (const double*)vIn;
  AddXY( ptr[0], ptr[1] );
}

// DataSet_Mesh::WriteBuffer()
void DataSet_Mesh::WriteBuffer(CpptrajFile &cbuffer, SizeArray const& pIn) const {
  if (pIn[0] >= mesh_x_.size())
    cbuffer.Printf(format_.fmt(), 0.0);
  else
    cbuffer.Printf(format_.fmt(), mesh_y_[pIn[0]]);
}
#ifdef MPI
// DataSet_double::Sync()
int DataSet_Mesh::Sync(size_t total, std::vector<int> const& rank_frames,
                       Parallel::Comm const& commIn)
{
  if (commIn.Size()==1) return 0;
  if (commIn.Master()) {
    // Resize for total number of frames.
    mesh_x_.resize( total );
    mesh_y_.resize( total );
    double* x_endptr = &(mesh_x_[0]) + rank_frames[0];
    double* y_endptr = &(mesh_y_[0]) + rank_frames[0];
    // Receive data from each rank
    for (int rank = 1; rank < commIn.Size(); rank++) {
      commIn.SendMaster( x_endptr, rank_frames[rank], rank, MPI_DOUBLE );
      commIn.SendMaster( y_endptr, rank_frames[rank], rank, MPI_DOUBLE );
      x_endptr += rank_frames[rank];
      y_endptr += rank_frames[rank];
    }
  } else {// Send data to master //TODO adjust for repeated additions?
    commIn.SendMaster( &(mesh_x_[0]), mesh_x_.size(), commIn.Rank(), MPI_DOUBLE );
    commIn.SendMaster( &(mesh_y_[0]), mesh_y_.size(), commIn.Rank(), MPI_DOUBLE );
  }
  return 0;
}
#endif

// -----------------------------------------------------------------------------
int DataSet_Mesh::Append(DataSet* dsIn) {
  if (dsIn->Empty()) return 0;
  if (dsIn->Group() != SCALAR_1D) return 1;
  if (dsIn->Type() == XYMESH) {
    std::vector<double> const& xIn = ((DataSet_Mesh*)dsIn)->mesh_x_;
    std::vector<double> const& yIn = ((DataSet_Mesh*)dsIn)->mesh_y_;
    size_t oldsize = Size();
    mesh_x_.resize( oldsize + xIn.size() );
    mesh_y_.resize( oldsize + yIn.size() );
    std::copy( xIn.begin(), xIn.end(), mesh_x_.begin() + oldsize );
    std::copy( yIn.begin(), yIn.end(), mesh_y_.begin() + oldsize );
  } else {
    DataSet_1D const& ds = static_cast<DataSet_1D const&>( *dsIn );
    for (unsigned int i = 0; i != ds.Size(); i++)
      AddXY( ds.Xcrd(i), ds.Dval(i) );
  }
  return 0;
}

// DataSet_Mesh::CalculateMeshX()
void DataSet_Mesh::CalculateMeshX(int sizeIn, double ti, double tf) {
  mesh_x_.resize( sizeIn, 0 );
  mesh_y_.resize( sizeIn, 0 );
  double s = (ti + tf)/2;
  double d = (tf - ti)/2;
  for (int i = 0; i < sizeIn; i++)
    mesh_x_[i] = s + d*((double) (2*i + 1 - sizeIn)/(sizeIn - 1));
  // Update dimension
  SetDim(Dimension::X, Dimension( ti, (tf - ti) / (double)(sizeIn - 1), Dim(0).Label() ));
}

// DataSet_Mesh::SetMeshXY()
int DataSet_Mesh::SetMeshXY(DataSet_1D const& dsIn) {
  mesh_x_.resize( dsIn.Size() );
  mesh_y_.resize( dsIn.Size() );
  for (unsigned int i = 0; i < dsIn.Size(); i++) {
    mesh_x_[i] = dsIn.Xcrd(i);
    mesh_y_[i] = dsIn.Dval(i);
  }
  // Update dimension
  SetDim(Dimension::X, Dimension( dsIn.Dim(0).Min(), dsIn.Dim(0).Step(), Dim(0).Label() ));
  return 0;
}

// ---------- Cubic Spline Routines --------------------------------------------
// DataSet_Mesh::SetSplinedMeshY()
/** Assumes mesh X values already set with CalculateMeshX. */
int DataSet_Mesh::SetSplinedMeshY(std::vector<double> const& x, std::vector<double> const& y) {
  if (x.size() != y.size()) {
    mprinterr("Error: X size (%zu) != Y size (%zu)\n", x.size(), y.size());
    return 1;
  }
  // No point if 1 or less values
  if (x.size() < 2) {
    mprinterr("Error: Requires > 1 values (%zu specified).\n", x.size());
    return 1;
  }
  cspline_.CubicSpline_Coeff(x, y);
  mesh_y_ = cspline_.CubicSpline_Eval(x, y, mesh_x_);
  return 0;
}

// DataSet_Mesh::SetSplinedMesh()
/** Assumes mesh X values already set with CalculateMeshX. */
int DataSet_Mesh::SetSplinedMesh(DataSet_1D const& dsIn)
{
  if (dsIn.Size() < 2) {
    mprinterr("Error: Requires > 1 values (%zu specified).\n", dsIn.Size());
    return 1;
  }
  // Create X and Y values for dsIn
  std::vector<double> x, y;
  x.reserve( dsIn.Size() );
  y.reserve( dsIn.Size() );
  for (int i = 0; i < (int)dsIn.Size(); i++) {
    x.push_back( dsIn.Xcrd( i ) );
    y.push_back( dsIn.Dval( i ) );
  }
  cspline_.CubicSpline_Coeff(x, y);
  mesh_y_ = cspline_.CubicSpline_Eval(x, y, mesh_x_);
  return 0;
}

// ---------- Regression -------------------------------------------------------
// DataSet_Mesh::SingleExponentialRegression()
int DataSet_Mesh::SingleExpRegression(double& slope, double& intercept,
                                      double& correl, CpptrajFile* out )
{
  std::vector<double> yorig = mesh_y_;
  for (unsigned int i = 0; i != mesh_y_.size(); i++)
  {
    if (mesh_y_[i] <= 0.0) {
      mprinterr("Error: '%s' Cannot perform exp. regression; set contains value <= 0\n", legend());
      mesh_y_ = yorig;
      return 1;
    }
    mesh_y_[i] = log( mesh_y_[i] );
  }
  int err = LinearRegression(slope, intercept, correl, out);
  // Restore original Y values
  mesh_y_ = yorig;
  return err;
}
