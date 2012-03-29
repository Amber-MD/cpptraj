#include "ParmIO.h"

// CONSTRUCTOR
ParmIO::ParmIO() { }

// Copy Constructor
ParmIO::ParmIO(const ParmIO &rhs) :
  CpptrajFile(rhs)
{ }

// Assignment
ParmIO &ParmIO::operator=(const ParmIO &rhs) {
  // Self
  if (this == &rhs) return *this;
  // Base Class
  CpptrajFile::operator=(rhs);
  // Deallocate
  // Allocate and copy
  return *this;
}

// ParmIO::SetDebug
void ParmIO::SetDebug(int debugIn) {
  debug_ = debugIn;
}

