#include "DataSet_LeapOpts.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
DataSet_LeapOpts::DataSet_LeapOpts() :
  // 0 dim indicates DataSet-specific write
  DataSet(LEAPOPTS, GENERIC, TextFormat(TextFormat::STRING, 12, 0), 0),
  pbradii_(Cpptraj::Parm::MBONDI),
  scee_(1.2), // AMBER default
  scnb_(2.0), // AMBER default
  dipoleDampFactor_(0),
  ipol_(0),
  flexibleWater_(false)
{

}

/** Set default GB radii from keyword. */
int DataSet_LeapOpts::SetGbRadii(std::string const& keyword) {
  Cpptraj::Parm::GB_RadiiType radType = Cpptraj::Parm::GbTypeFromKey( keyword );
  if (radType == Cpptraj::Parm::UNKNOWN_GB) {
    mprinterr("Error: Unrecognized GB radii type: %s\n", keyword.c_str());
    return 1;
  }
  pbradii_ = radType;
  mprintf("\tSet default GB radii: %s\n", Cpptraj::Parm::GbAmberFlag(pbradii_).c_str());
  return 0;
}

/// Set double option, warn if negative
static inline int set_dopt(const char* desc, double& opt, double newVal)
{
  if (newVal < 0.0) {
    mprintf("Warning: Illegal value: %s can not be < 0.0;\n"
            "Warning: keeping value of %g\n", desc, opt);
    return 1;
  } else {
    mprintf("\tSetting %s to %g\n", desc, newVal);
    opt = newVal;
  }
  return 0;
}

/** Set default SCEE */
int DataSet_LeapOpts::SetSCEE(double sceeIn) { return set_dopt("1-4 SCEE factor", scee_, sceeIn); }

/** Set default SCNB */
int DataSet_LeapOpts::SetSCNB(double scnbIn) { return set_dopt("1-4 SCNB factor", scnb_, scnbIn); }

/** Set dipole damping factor */
int DataSet_LeapOpts::SetDipoleDampFactor(double ddfIn) { return set_dopt("Dipole damping factor", dipoleDampFactor_, ddfIn); }

/** Set IPOL */
int DataSet_LeapOpts::SetIpol(int ipolIn) {
  if (ipolIn < 0 || ipolIn > 4) {
    mprintf("Warning: Only IPOL = 0 to 4 is supported, resetting IPOL to %i.\n", ipol_);
    return 1;
  }
  if (ipol_ > 0) {
    mprintf("Warning: IPOL has already been set to %i by a previous force field.\n", ipol_);
    return 1;
  }
  mprintf("\tSetting IPOL to %i\n", ipolIn);
  ipol_ = ipolIn;
  return 0;
}

/** Set flexible water */
int DataSet_LeapOpts::SetFlexibleWater(bool fIn) {
  if (fIn)
    mprintf("\tSetting flexible water to ON\n");
  else
    mprintf("\tSetting flexible water to OFF\n");
  flexibleWater_ = fIn;
  return 0;
}
