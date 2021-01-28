#include "EnergyArray.h"
#include "CpptrajStdio.h"
#include "CpptrajFile.h"

/** CONSTRUCTOR */
EnergyArray::EnergyArray() :
  ene_((int)N_E_TERMS, 0.0)
{}

/** Energy term labels. Keep in sync with Type. */
const char* EnergyArray::TypeStr_[] = { "Bond", "Angle", "Dihedral", "VDW14", "Q14", "VDW", "Elec", "OpenMM", 0 };

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

/** Print active terms to the given file. */
void EnergyArray::PrintActiveTerms(CpptrajFile& outfile, bool includeSingleTerms)
const
{
  if (activeTerms_.empty()) return;
  if (!includeSingleTerms && activeTerms_.size() < 2) return;
  for (Tarray::const_iterator it = activeTerms_.begin(); it != activeTerms_.end(); ++it)
    outfile.Printf(" %12.4E", ene_[*it]);
}

/** Print labels of active terms to the given file. */
void EnergyArray::PrintActiveLabels(CpptrajFile& outfile, bool includeSingleTerms)
const
{
  if (activeTerms_.empty()) return;
  if (!includeSingleTerms && activeTerms_.size() < 2) return;
  for (Tarray::const_iterator it = activeTerms_.begin(); it != activeTerms_.end(); ++it)
    outfile.Printf(" %12s", TypeStr_[*it]);
}

/** Print active term labels to STDOUT. */
void EnergyArray::PrintActiveTerms() const {
  for (Tarray::const_iterator it = activeTerms_.begin(); it != activeTerms_.end(); ++it)
    mprintf(" %s", TypeStr_[*it]);
}
