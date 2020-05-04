#include "EnergyArray.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
EnergyArray::EnergyArray() :
  ene_((int)N_E_TERMS, 0.0)
{}

/** Energy term labels. Keep in sync with Type. */
const char* EnergyArray::TypeStr_[] = { "Bond", "VDW", "Elec", 0 };

const char* EnergyArray::label(Type t) const {
  return TypeStr_[(int)t];
}

/** Specify that an energy term will be calculated.
  * \return Pointer to position in energy array of specified term.
  */
double* EnergyArray::AddType(Type typeIn) {
  for (Tarray::const_iterator it = activeTerms_.begin(); it != activeTerms_.end(); ++it)
  {
    if (*it == typeIn) {
      mprinterr("Error: Energy term %s already present.\n", TypeStr_[typeIn]);
      return 0;
    }
  }
  activeTerms_.push_back( typeIn );
  return ( (&ene_[0]) + ((int)typeIn) );
}
