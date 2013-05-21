#include "Array1D.h"
#include "CpptrajStdio.h"
// CONSTRUCTOR
Array1D::Array1D(DataSetList const& SetList) {
  AddDataSets( SetList );
  if (array_.empty())
    mprinterr("Internal Error: No 1D data sets present.");
}

// Array1D::push_back()
int Array1D::push_back( DataSet_1D* const& val ) {
  if (val->Ndim() == 1)
    array_.push_back( val );
  else
    return 1;
  return 0;
}

// Array1D::AddDataSets()
int Array1D::AddDataSets(DataSetList const& SetList) {
  for (DataSetList::const_iterator ds = SetList.begin(); ds != SetList.end(); ++ds)
    if ( push_back( (DataSet_1D*)*ds ) ) {
      mprinterr("Error: Only 1D data sets allowed.");
      array_.clear();
      return 1;
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
  // Sort input datasets
  input_dsl.sort();
  // Add to main list
  array_.clear();
  if (AddDataSets( input_dsl ))
    return 1;
  return 0;
}

// Array1D::DetermineMax() 
size_t Array1D::DetermineMax() const {
  size_t maxFrames = 0L;
  for (std::vector<DataSet_1D*>::const_iterator set = array_.begin(); set != array_.end(); ++set)
    if ( (*set)->Size() > maxFrames )
      maxFrames = (*set)->Size();
  return maxFrames;
}
