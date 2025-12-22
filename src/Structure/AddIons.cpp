#include "AddIons.h"
#include "../CpptrajStdio.h"

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
AddIons::AddIons() :
  Nion1_(-1),
  Nion2_(-1),
  debug_(0),
  separation_(0.0)
{}

/** Initialize */
int AddIons::InitAddIons(std::string const& unitNameIn,
                         std::string const& ion1nameIn, int Nion1,
                         std::string const& ion2nameIn, int Nion2,
                         double separationIn, int seedIn, int debugIn)
{
  debug_ = debugIn;

  if (unitNameIn.empty()) {
    mprinterr("Internal Error: AddIons::InitAddIons() called with empty unit name.\n");
    return 1;
  }
  unitname_ = unitNameIn;

  if (ion1nameIn.empty()) {
    mprinterr("Internal Error: AddIons::InitAddIons() called with empty ion name.\n");
    return 1;
  }
  ion1name_ = ion1nameIn;
  if (Nion1 < 0) {
    mprinterr("Error: Number of %s ions cannot be less than 0 (%i)\n", ion1nameIn.c_str(), Nion1);
    return 1;
  }
  Nion1_ = Nion1;

  if (ion2nameIn.empty()) {
    ion2name_.clear();
    Nion2_ = -1;
  } else {
    ion2name_ = ion2nameIn;
    if (Nion2 < 0) {
      mprinterr("Error: Number of %s ions cannot be less than 0 (%i)\n", ion2nameIn.c_str(), Nion2);
      return 1;
    }
  }

  separation_ = separationIn;
  RNG_.rn_set( seedIn );

  return 0;
}
