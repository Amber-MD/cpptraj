#include "DataSet_Coords_FRM.h"

/** CONSTRUCTOR */
DataSet_Coords_FRM::DataSet_Coords_FRM() {}

#ifdef MPI
// FIXME Sync
#endif

/** Add a single element. */
void DataSet_Coords_FRM::Add(size_t idx, const void* frmPtr) {
  if (idx > frames_.size())
    frames_.resize( idx );
  // Insert at end
  frames_.push_back( *((Frame*)frmPtr) );
}


