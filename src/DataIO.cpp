#include "DataIO.h"
#include "CpptrajStdio.h"
#include "DataSet_MatrixDbl.h"

// DataIO::CheckValidFor()
bool DataIO::CheckValidFor( DataSet const& dataIn ) const {
  if (valid1d_ && dataIn.Ndim() == 1) return true; 
  if (valid2d_ && dataIn.Ndim() == 2) return true; 
  if (valid3d_ && dataIn.Ndim() == 3) return true;
  for (std::vector<DataSet::DataType>::const_iterator dt = valid_.begin(); 
                                                      dt != valid_.end(); ++dt)
    if (dataIn.Type() == *dt) return true;
  return false;
}

int DataIO::CheckAllDims(DataSetList const& array, unsigned int tgtDim) {
  for (DataSetList::const_iterator set = array.begin(); set != array.end(); ++set)
  {
    if ( (*set)->Ndim() != tgtDim ) {
      mprinterr("Error: Set '%s' dimension is %i, expected only %iD.\n",
                (*set)->legend(), (*set)->Ndim(), tgtDim);
      return 1;
    }
  }
  return 0;
}

int DataIO::CheckXDimension(DataSetList const& array) {
  if (array.empty()) return 0; // FIXME return error?
  int err = 0;
  Dimension const& Xdim = static_cast<Dimension const&>(array[0]->Dim(0));
  for (DataSetList::const_iterator set = array.begin(); set != array.end(); ++set)
  {
    if ((*set)->Dim(0) != Xdim) {
      mprinterr("Error: X Dimension of %s != %s\n", (*set)->legend(),
                array[0]->legend());
      mprinterr("Error:  %s: Min=%f Step=%f\n", (*set)->legend(),
                (*set)->Dim(0).Min(), (*set)->Dim(0).Step());
      mprinterr("Error:  %s: Min=%f Step=%f\n", array[0]->legend(),
                Xdim.Min(), Xdim.Step());
      ++err;
    }
  }
  return err;
}

size_t DataIO::DetermineMax(DataSetList const& array) {
  size_t maxSize = 0L;
  for (DataSetList::const_iterator set = array.begin(); set != array.end(); ++set)
    if ( (*set)->Size() > maxSize )
      maxSize = (*set)->Size();
  return maxSize;
}

/** Given a flattened matrix (row-major order) with given # rows and columns,
  * determine if it is symmetric and allocate into the given DataSetList
  * accordingly.
  */
DataSet* DataIO::DetermineMatrixType(std::vector<double> const& matrixArray, int nrows, int ncols,
                                DataSetList& DSL, std::string const& dsname)
{
  DataSet* ds = DSL.AddSet(DataSet::MATRIX_DBL, dsname, "Mat");
  if (ds == 0) return 0;
  DataSet_MatrixDbl& Mat = static_cast<DataSet_MatrixDbl&>( *ds );
  //ds->SetupMeta().SetScalarType( MetaData::DIST ); // TODO: FIXME Allow type keywords
  bool isSymmetric = false;
  if (ncols == nrows) {
    isSymmetric = true;
    // Check if matrix is symmetric
    for (int row = 0; row < nrows; row++) {
      for (int col = row + 1; col < ncols; col++) {
        if ( matrixArray[ (row * ncols) + col ] != matrixArray[ (col * ncols) + row ] ) {
          isSymmetric = false;
          break;
        }
      }
      if (!isSymmetric) break;
    }
  }
  if (isSymmetric) {
    mprintf("\tSymmetric matrix detected.\n");
    if (Mat.AllocateHalf(ncols)) {
      mprinterr("Error: Could not allocate memory for set '%s'\n", ds->legend());
      DSL.RemoveSet( ds );
      return 0;
    }
    for (int row = 0; row < nrows; row++)
      for (int col = row; col < ncols; col++)
        Mat.AddElement( matrixArray[ (row * ncols) + col ] );
  } else {
    DataSet::SizeArray dims(2);
    dims[0] = ncols;
    dims[1] = nrows;
    ds->Allocate( dims );
    std::copy( matrixArray.begin(), matrixArray.end(), Mat.begin() );
  }
  return ds;
}
