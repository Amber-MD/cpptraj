#include "DataIO.h"
#include "Constants.h"
#include "CpptrajStdio.h"
#include "DataSet_1D.h"
#include "DataSet_MatrixDbl.h"
#include <algorithm> // std::min

// CONSTRUCTOR
DataIO::DataIO() :
  debug_(0),
  xcol_fmt_(TextFormat::DOUBLE), // default
  xcol_width_(8),                // default
  xcol_prec_(3),                 // default
  x_prec_set_(false),
  valid1d_(false),
  valid2d_(false),
  valid3d_(false)
{}

/** CONSTRUCTOR - Set valid for 1d, 2d, and/or 3d data sets */
DataIO::DataIO(bool v1, bool v2, bool v3) :
  debug_(0),
  xcol_fmt_(TextFormat::DOUBLE), // default
  xcol_width_(8),                // default
  xcol_prec_(3),                 // default
  x_prec_set_(false),
  valid1d_(v1),
  valid2d_(v2),
  valid3d_(v3)
{}

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

/** \return 1 if any set in the given array does not have the given number of dimensions, 0 otherwise.
  */
int DataIO::CheckAllDims(DataSetList const& array, unsigned int tgtDim) {
  for (DataSetList::const_iterator set = array.begin(); set != array.end(); ++set)
  {
    if ( (*set)->Ndim() != tgtDim ) {
      mprinterr("Error: Set '%s' dimension is %zu, expected only %uD.\n",
                (*set)->legend(), (*set)->Ndim(), tgtDim);
      return 1;
    }
  }
  return 0;
}
/*
/// \return true if floating point values are equivalent. Used to check X values.
static inline bool FEQ(double v1, double v2) {
  double delta = v1 - v2;
  if (delta < 0.0) delta = -delta;
  return (delta < Constants::SMALL);
}

/// \return true if floating point values are not equivalent. Used to check X values.
static inline bool FNE(double v1, double v2) {
  double delta = v1 - v2;
  if (delta < 0.0) delta = -delta;
  return (delta > Constants::SMALL);
}*/

/// \return X coordinate of given set.
static inline double xCoordVal(DataSet const& set, unsigned int idx) {
  if (set.Group() == DataSet::SCALAR_1D) {
    DataSet_1D const& s1 = static_cast<DataSet_1D const&>( set );
    return s1.Xcrd(idx);
  } else
    return set.Dim(0).Coord(idx);
}

/** Check if X dimension of 2 given sets matches. */
bool DataIO::xDimMatch(DataSet const& ref, DataSet const& set) {
  if (ref.Ndim() < 1 || set.Ndim() < 1) return false;
  // First check if the reference set is the index set for incoming set.
  if ( set.DimIndexSet(0) != 0 && set.DimIndexSet(0) == (DataSet const*)(&ref) ) {
    //mprintf("DEBUG: %s is the X index set for %s\n", ref.legend(), set.legend());
    return true;
  } else if (ref.Type() == DataSet::XYMESH || set.Type() == DataSet::XYMESH) {
    // Either one of both sets is a mesh. Need to explicitly check X values.
    if (ref.Size() == 0 && set.Size() == 0)
      // No values in either; automatic match
      return true;
    // Start from the end since that is most likely to not match
    unsigned int endval = std::min(set.Size(), ref.Size());
    if (endval == 0)
      // One of the sets is empty; automatic non-match TODO should just match?
      return false;
    //for (int idx = endval - 1; idx >= 0; idx--)
    for (unsigned int idx = 0; idx < endval; idx++)
    {
      double refXval = xCoordVal(ref, idx);
      double setXval = xCoordVal(set, idx);
      if (FNE( refXval, setXval )) {
        mprintf("Warning: X coordinate of %s (%f) != %s (%f) at index %u\n",
                set.legend(), setXval, ref.legend(), refXval, idx);
        return false;
      }
    }
  } else if (set.Dim(0) != ref.Dim(0)) {
    mprintf("Warning: X Dimension of %s (Min=%f Step=%f) != %s (Min=%f Step=%f)\n",
            set.legend(), set.Dim(0).Min(), set.Dim(0).Step(),
            ref.legend(), ref.Dim(0).Min(), ref.Dim(0).Step());
    return false;
  }
  return true;
}

/** Ensure that the X dimension of each set in the array matches.
  * \return An array of DataSetLists where the sets in each
  *         DataSetList all have matching X dimensions.
  */
DataIO::DSLarray DataIO::CheckXDimension(DataSetList const& array) {
  DSLarray OUT;
  if (array.empty()) return OUT; // FIXME return error?
  //Xmatch.clear();
  //Xmatch.resize( array.size(), -1 );
  std::vector<int> Xmatch(array.size(), -1);

  // Loop over all sets
  unsigned int NsetsChecked = 0;
//  int err = 0;
//  Dimension const& Xdim = static_cast<Dimension const&>(array[0]->Dim(0));
  while (NsetsChecked < array.size()) {
    // Find the first unmatched set
    unsigned int currentIndex = 0; 
    while (currentIndex < array.size() && Xmatch[currentIndex] != -1) ++currentIndex;
    if (currentIndex >= array.size()) break;
    DataSet const& currentRef = *(array[currentIndex]);
//    mprintf("DEBUG: Current ref: %s\n", currentRef.legend());
    Xmatch[currentIndex] = currentIndex;
    OUT.push_back( DataSetList() );
    OUT.back().AddCopyOfSet( array[currentIndex] );
    // Loop over remaining sets
//    mprintf("DEBUG: Starting at set %u\n", currentIndex+1);
    for (unsigned int setIdx = currentIndex+1; setIdx < array.size(); setIdx++)
    {
      DataSet const& currentSet = *array[setIdx];
      if (Xmatch[setIdx] == -1) {
//        mprintf("DEBUG:\t\tChecking %s", currentSet.legend());
        // Check if currentSet Xdim matches the current reference
        if (xDimMatch( currentRef, currentSet )) {
//          mprintf(" MATCH!");
          Xmatch[setIdx] = (int)currentIndex;
          OUT.back().AddCopyOfSet( array[setIdx] );
          NsetsChecked++;
        }
//        mprintf("\n");
      }
    }
  }
  // DEBUG
/*  mprintf("DEBUG: Xdim match indices:\n");
  for (unsigned int idx = 0; idx != array.size(); idx++)
    mprintf("DEBUG:\t%20s %i\n", array[idx]->legend(), Xmatch[idx]);
  mprintf("DEBUG: Sets grouped by Xdim:\n");
  for (DSLarray::const_iterator it = OUT.begin(); it != OUT.end(); ++it) {
    mprintf("\t");
    for (DataSetList::const_iterator jt = it->begin(); jt != it->end(); ++jt)
      mprintf(" %s", (*jt)->legend());
    mprintf("\n");
  }*/
  if (OUT.size() > 1) {
    mprintf("Warning: Not all sets had matching X dimensions. They will be grouped in the file like so:\n");
    for (DSLarray::const_iterator it = OUT.begin(); it != OUT.end(); ++it) {
      mprintf("\t");
      for (DataSetList::const_iterator jt = it->begin(); jt != it->end(); ++jt)
        mprintf(" %s", (*jt)->legend());
      mprintf("\n");
    }
  }
//  DataSetList::const_iterator set = array.begin();
//  for (DataSetList::const_iterator set = array.begin(); set != array.end(); ++set)
//  {
//    Dimension const& currentRefXdim = static_cast<Dimension const&>( array[currentMatchIndex]->Dim(0) );
//  }
//  return err;
  return OUT;
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
