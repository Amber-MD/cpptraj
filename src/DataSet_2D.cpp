#include "DataSet_2D.h"

/// Strings corresponding to MatrixType enumerated type.
const char* DataSet_2D::MatrixTypeString[] = {
  "UNDEFINED",
  "distance matrix",
  "covar matrix",
  "mass weighted covar matrix",
  "correlation matrix",
  "distance covar matrix",
  "idea matrix",
  "ired matrix"
};

const char* DataSet_2D::MatrixOutputString[] = {
  "UNKNOWN",
  "DIST",
  "COVAR",
  "MWCOVAR",
  "CORREL",
  "DISTCOVAR",
  "IDEA",
  "IRED"
};
