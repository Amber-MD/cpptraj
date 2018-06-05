#include "DataSet_Parameters.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
DataSet_Parameters::DataSet_Parameters() :
  DataSet(PARAMETERS, GENERIC, TextFormat(), 0)
{}

size_t DataSet_Parameters::Size() const {
  return (atomTypes_.Size() +
          bondParm_.size() +
          angleParm_.size() +
          ubParm_.size() +
          dihParm_.size() +
          impParm_.size());
}

void DataSet_Parameters::Info() const {
  if (Size() > 0) {
    mprintf(" (");
    if (atomTypes_.Size() > 0) mprintf(" types=%zu", atomTypes_.Size());
    if (bondParm_.size() > 0) mprintf(" bnd=%zu", bondParm_.size());
    if (angleParm_.size() > 0) mprintf(" ang=%zu", angleParm_.size());
    if (ubParm_.size() > 0) mprintf(" UB=%zu", ubParm_.size());
    if (dihParm_.size() > 0) mprintf(" dih=%zu", dihParm_.size());
    if (impParm_.size() > 0) mprintf(" imp=%zu", impParm_.size());
    mprintf(" )");
  }
}
