#include "Units.h"
#include "Constants.h"
#include "CpptrajStdio.h"
#include <cctype> // tolower

/** \return Unit type from name. */
Cpptraj::Units::Type Cpptraj::Units::TypeFromName(std::string const& nameIn) {
  std::string name;
  for (std::string::const_iterator it = nameIn.begin(); it != nameIn.end(); ++it)
    name += tolower( *it );
  if (name == "angstroms"  || name == "angstrom"  || name == "ang") return ANG;
  if (name == "nanometers" || name == "nanometer" || name == "nm") return NM;
  mprintf("Warning: Unrecognized units: %s\n", nameIn.c_str());
  return UNKNOWN_UNITS;
}

/** Set the conversion factor required to convert the first unit
  * into the second unit via multiplication.
  */
int Cpptraj::Units::SetConversionFactor(double& fac,
                                        std::string const& fromName, std::string const& toName)
{
  fac = 0;
  Type fromUnits = TypeFromName(fromName);
  if (fromUnits == UNKNOWN_UNITS) return 1;
  Type toUnits = TypeFromName(toName);
  if (toUnits == UNKNOWN_UNITS) return 1;

  // If units match no conversion factor needed
  if (fromUnits == toUnits) {
    fac = 1.0;
    return 0;
  }

  if (fromUnits == ANG) {
    // from angstroms to X
    if (toUnits == NM) {
      fac = Constants::ANG_TO_NM;
      return 0;
    }
  } else if (fromUnits == NM) {
    // from nanometers to X
    if (toUnits == ANG) {
      fac = Constants::NM_TO_ANG;
      return 0;
    }
  }
  mprintf("Warning: No conversion from '%s' to '%s'\n", fromName.c_str(), toName.c_str());
  return 1;
}
