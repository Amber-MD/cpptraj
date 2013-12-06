#ifndef INC_CONSTANTS_H
#define INC_CONSTANTS_H
/*! \file Constants.h
    \brief Useful Physical constants

    The constants here are stored in as high precision as possible for
    double. Certain constants may have larger than 52 bit fraction,
    which may result in some small precision loss.
 */
namespace Constants {
  // Various incarnations of PI
  const double PI           = 3.1415926535897932384626433832795; 
  const double TWOPI        = 6.2831853071795864769252867665590; 
  const double FOURPI       = 12.566370614359172953850573533118; 
  const double FOURTHIRDSPI = 4.1887902047863909846168578443727; 
  const double FOURFIFTHSPI = 2.5132741228718345907701147066236; 
  const double PIOVER2      = 1.5707963267948966192313216916398; 
  // Convert degrees <-> radians
  const double DEGRAD       = 0.017453292519943295769236907684886; 
  const double RADDEG       =   57.295779513082320876798154814105; 
  // For checking floating point zero
  const double SMALL        = 0.00000000000001;
  /// Gas constant in J/mol*K
  const double GASK_J       = 8.3144621;
  /// Gas constant in kcal/mol*K
  const double GASK_KCAL    = 0.0019872041;
  // Convert electron charge <-> Amber units (w/ prefactor)
  const double ELECTOAMBER  = 18.2223;
  const double AMBERTOELEC  = 1.0 / ELECTOAMBER;
}
#endif
