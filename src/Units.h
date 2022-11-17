#ifndef INC_UNITS_H
#define INC_UNITS_H
#include <string>
namespace Cpptraj {
namespace Units {

enum Type { ANG = 0, ///< Angstroms
            NM,      ///< nanometers
            UNKNOWN_UNITS
           };

/// \return Unit type from name.
Type TypeFromName(std::string const&);

/// Set multiplicative conversion factor from first unit to second unit
int SetConversionFactor(double&, std::string const&, std::string const&);

}
}
#endif
