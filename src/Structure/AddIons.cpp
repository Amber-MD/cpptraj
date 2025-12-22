#include "AddIons.h"
#include "../CpptrajStdio.h"
#include "../Topology.h"

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
AddIons::AddIons() :
  Nion1_(-1),
  Nion2_(-1),
  debug_(0),
  separation_(0.0)
{}

/** Initialize */
int AddIons::InitAddIons(std::string const& ion1nameIn, int Nion1,
                         std::string const& ion2nameIn, int Nion2,
                         double separationIn, int seedIn, int debugIn)
{
  debug_ = debugIn;

  if (ion1nameIn.empty()) {
    mprinterr("Error: Must specify at least one ion name.\n");
    return 1;
  }
  ion1name_ = ion1nameIn;
  if (Nion1 < 0) {
    //mprinterr("Error: Number of %s ions cannot be less than 0 (%i)\n", ion1nameIn.c_str(), Nion1);
    //return 1;
    Nion1_ = 0;
  } else
    Nion1_ = Nion1;

  if (ion2nameIn.empty()) {
    ion2name_.clear();
    Nion2_ = -1;
  } else {
    ion2name_ = ion2nameIn;
    if (Nion2 < 1) {
      mprinterr("Error: Number of second %s ions cannot be less than 1 (%i)\n", ion2nameIn.c_str(), Nion2);
      return 1;
    }
  }

  separation_ = separationIn;
  if (separation_ < 0.0) {
    mprinterr("Error: Separation must be >= 0.0 (%f)\n", separation_);
    return 1;
  }

  RNG_.rn_set( seedIn );

  return 0;
}

/** Print setup info to stdout */
void AddIons::PrintAddIonsInfo() const {
  if (ion1name_.empty()) return;
  if (Nion1_ < 1)
    mprintf("\tAdding enough %s ions to neutralize.\n", ion1name_.c_str());
  else
    mprintf("\tAdding %i %s ions.\n", Nion1_, ion1name_.c_str());
  if (!ion2name_.empty())
    mprintf("\tAdding %i %s ions.\n", Nion2_, ion2name_.c_str());
  mprintf("\tMinimum ion separation is %g Ang.\n", separation_);
  mprintf("\tIon RNG seed: %i\n", RNG_.Seed());
}

/** Add ions randomly, replacing solvent molecules. */
int AddIons::AddIonsRand(Topology& topOut, Frame& frameOut, DataSetList const& DSL,
                         Cpptraj::Parm::ParameterSet const& set0)
const
{
  // Unit total charge
  double totalCharge = topOut.TotalCharge();
  if ( fabs(totalCharge) < Constants::SMALL ) {
    mprintf("Warning: %s has a total charge of 0.\n", topOut.c_str());
    if (Nion1_ < 1) {
      mprintf("Warning: Cannot neutralize.\n");
      return 0;
    }
    mprintf("Warning: Adding the ions anyway.\n");
  } else {
    mprintf("\tTotal charge on %s is %g\n", topOut.c_str(), totalCharge);
  }

  return 0;
}
