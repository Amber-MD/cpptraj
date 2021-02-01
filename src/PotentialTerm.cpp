#include "PotentialTerm.h"

const char* PotentialTerm::typeStr_[] = {
  "Bond",
  "Angle",
  "Dihedral",
  "SimpleNonbond",
  "OpenMM",
  0
};

const char* PotentialTerm::TypeStr(Type typeIn) {
  return typeStr_[typeIn];
}
