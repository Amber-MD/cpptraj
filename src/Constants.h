#ifndef INC_CONSTANTS_H
#define INC_CONSTANTS_H
/*! \file Constants.h
    \brief Useful Physical constants

    The constants here are stored in as high precision as possible for
    double. Certain constants may have larger than 52 bit fraction,
    which may result in some small precision loss. The number of digits
    chosen in those cases (all PI-related) was done so using the criterion
    that the conversion back to PI matches within roundoff (e.g.
    TWOPI / 2 = 3.141592653589793).
 */
namespace Constants {
  // Various incarnations of PI
  const double PI           = 3.141592653589793; 
  const double TWOPI        = 6.283185307179586;  // may have precision loss
  const double FOURPI       = 12.566370614359172; // may have precision loss
  const double FOURTHIRDSPI = 4.1887902047863909; // may have precision loss 
  const double FOURFIFTHSPI = 2.5132741228718345; // may have precision loss
  const double PIOVER2      = 1.5707963267948966; // may have precision loss 
  /// Convert degrees to radians
  const double DEGRAD       = 0.017453292519943295; // may have precision loss 
  /// Convert radians to degrees
  const double RADDEG       =   57.29577951308232;  // may have precision loss
  /// For checking floating point zero
  const double SMALL        = 0.00000000000001;
  /// Gas constant in J/mol*K
  const double GASK_J       = 8.3144621;
  /// Gas constant in kcal/mol*K
  const double GASK_KCAL    = 0.0019872041;
  // Avogadro constant
  const double NA = 6.02214129e23;
  // Speed of light (m/s)
  const double C0 = 299792458;
  /// Convert atomic mass unit (amu) to kg
  const double AMU_TO_KG = 1.660539e-27;
  /// Convert Angstroms to nanometers
  const double ANG_TO_NM = 0.1;
  /// Convert nanometers to Angstroms
  const double NM_TO_ANG = 10.0;
  /// Convert calories to Joules
  const double CAL_TO_J = 4.184;
  /// Convert Joules to calories
  const double J_TO_CAL = 1.0 / CAL_TO_J;
  /// Convert electron charge to Amber units (w/ prefactor)
  const double ELECTOAMBER  = 18.2223;
  /// Convert Amber charge to electron charge
  const double AMBERTOELEC  = 1.0 / ELECTOAMBER;
  /// Convert from Amber internal units of time (1/20.455 ps) to ps.
  /** Amber operates in kcal/mol units for energy, amu for masses,
    * and Angstoms for distances. For convenience when calculating KE from
    * velocity, the velocities have a conversion factor built in; as a result
    * the Amber unit of time is (1/20.455) ps. So e.g. to convert Amber
    * velocities from internal units to Ang/ps multiply by 20.455. The number
    * itself is derived from sqrt(1 / ((AMU_TO_KG * NA) / (1000 * CAL_TO_J))).
    */
  const double AMBERTIME_TO_PS = 20.455;
}
#endif
