#include "DataSet_Parameters.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
DataSet_Parameters::DataSet_Parameters() {}

void DataSet_Parameters::Debug() const {
  mprintf("Atom Types:\n");
  for (AtomTypeArray::const_iterator at = atomTypes_.begin(); at != atomTypes_.end(); ++at)
    mprintf("\t%s %i\n", *(at->first), at->second);
}
  
