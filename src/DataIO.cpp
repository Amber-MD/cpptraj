#include "DataIO.h"
#include "CpptrajStdio.h"

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
