#include <algorithm> // sort
#include "Array1D.h"
#include "CpptrajStdio.h"

// COPY CONSTRUCTOR
Array1D::Array1D(const Array1D& rhs) : array_(rhs.array_) {}

// ASSIGNMENT
Array1D& Array1D::operator=(const Array1D& rhs) {
  if (this == &rhs) return *this;
  array_ = rhs.array_;
  return *this;
}

// CONSTRUCTOR
Array1D::Array1D(DataSetList const& SetList) {
  AddDataSets( SetList );
  if (array_.empty())
    mprinterr("Internal Error: No 1D data sets present.\n");
}

// Array1D::push_back()
int Array1D::push_back( DataSet* val ) {
  if (val == 0) {
    mprinterr("Internal Error: Blank pointer passed to Array1D.\n");
    return 1;
  } else if (val->Group() != DataSet::SCALAR_1D) {
    mprintf("Warning: Cannot add '%s'; only 1D data sets allowed.\n", val->legend());
    return 0;
  }
  array_.push_back( (DataSet_1D*)val );
  return 0;
}

void Array1D::SortArray1D() {
  std::sort( array_.begin(), array_.end(), DataSet::DS_PtrCmp() );
}

// Array1D::AddDataSets()
int Array1D::AddDataSets(DataSetList const& SetList) {
  for (DataSetList::const_iterator ds = SetList.begin(); ds != SetList.end(); ++ds)
    if ( push_back( *ds ) ) {
      array_.clear();
      return 1;
    }
  return 0;
}

// Array1D::AddTorsionSets()
int Array1D::AddTorsionSets(DataSetList const& SetList) {
  // Ensure data sets are 1D and periodic
  for (DataSetList::const_iterator ds = SetList.begin(); ds != SetList.end(); ++ds) {
    if ( (*ds)->Meta().IsTorsionArray()) {
      if ( push_back( *ds ) ) {
        array_.clear();
        return 1;
      }
    } else
        mprintf("Warning: Set '%s' is not periodic, skipping.\n", (*ds)->legend());
  }
  return 0;
}

// Array1D::AddSetsFromArgs()
int Array1D::AddSetsFromArgs(ArgList const& dsetArgs, DataSetList const& DSLin) {
  DataSetList input_dsl;
  for (ArgList::const_iterator dsa = dsetArgs.begin(); dsa != dsetArgs.end(); ++dsa)
    input_dsl += DSLin.GetMultipleSets( *dsa );
  if (input_dsl.empty()) {
    mprinterr("Error: No data sets selected.\n");
    return 1;
  }
  // Add to main list
  array_.clear();
  if (AddDataSets( input_dsl ))
    return 1;
  return 0;
}
