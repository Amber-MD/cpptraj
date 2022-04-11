#include "DataSet_Parameters.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
DataSet_Parameters::DataSet_Parameters() :
  DataSet(PARAMETERS, GENERIC, TextFormat(), 0)
{}

size_t DataSet_Parameters::Size() const {
  return (AT().size() +
          BP().size() +
          AP().size() +
          UB().size() +
          DP().size() +
          IP().size());
}

void DataSet_Parameters::Info() const {
  if (Size() > 0) {
    mprintf(" (");
    if (AT().size() > 0) mprintf(" types=%zu", AT().size());
    if (NB().size() > 0) mprintf(" LJ pairs=%zu", NB().size());
    if (BP().size() > 0) mprintf(" bnd=%zu", BP().size());
    if (AP().size() > 0) mprintf(" ang=%zu", AP().size());
    if (UB().size() > 0) mprintf(" UB=%zu", UB().size());
    if (DP().size() > 0) mprintf(" dih=%zu", DP().size());
    if (IP().size() > 0) mprintf(" imp=%zu", IP().size());
    mprintf(" )");
  }
}
