#include "AddIons.h"
#include "../Parm/ParameterSet.h"
#include "../Parm/ParmHolder.h"
#include "../AtomType.h"
#include "../Constants.h"
#include "../CpptrajStdio.h"
#include "../DataSet_Coords.h"
#include "../DataSetList.h"
#include "../Topology.h"
#include <cmath> //fabs, lrint

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

/** Get ion unit box from DataSetList. */ // TODO consolidate with Solvate::GetSolventUnit()?
DataSet_Coords* AddIons::GetIonUnit(std::string const& ionname, DataSetList const& DSL) const {
  if (ionname.empty()) {
    mprinterr("Internal Error: AddIons::GetIonUnit() called before ion name set.\n");
    return 0;
  }
  DataSetList sets = DSL.SelectGroupSets( "*", DataSet::COORDINATES ); // TODO specific set type for units?
  // First try to match aspect, then match name
  DataSet_Coords* ionUnit = 0;
  // Aspect
  for (DataSetList::const_iterator it = sets.begin(); it != sets.end(); ++it)
  {
    DataSet_Coords* ds = static_cast<DataSet_Coords*>( *it );
    if (!ds->Meta().Aspect().empty()) {
      if (ionname == ds->Meta().Aspect()) {
        ionUnit = ds;
        break;
      }
    }
  }
  // Name
  if (ionUnit == 0) {
    for (DataSetList::const_iterator it = sets.begin(); it != sets.end(); ++it)
    {
      DataSet_Coords* ds = static_cast<DataSet_Coords*>( *it );
      if (ionname == ds->Meta().Name()) {
        ionUnit = ds;
        break;
      }
    }
  }
  if (ionUnit != 0)
    mprintf("\tIon unit: %s\n", ionUnit->legend());
  else
    mprinterr("Error: Could not get ion unit named %s\n", ionname.c_str());

  return ionUnit;
}

/// \return Array containing radii for every atom
static inline std::vector<double> GetAtomRadii(DataSet_Coords* crd, Cpptraj::Parm::ParmHolder<AtomType> const& newAtomTypeParams)
{
  static const double ATOM_DEFAULT_RADIUS = 1.5; // To match LEAP
  std::vector<double> OUT;
  OUT.reserve( crd->Top().Natom() );

  for (std::vector<Atom>::const_iterator at = crd->Top().begin(); at != crd->Top().end(); ++at)
  {
    bool found;
    TypeNameHolder atype( at->Type() );
    AtomType AT = newAtomTypeParams.FindParam( atype, found );
    if (found && AT.HasLJ())
      OUT.push_back( AT.LJ().Radius() );
    else {
      mprintf("Warning: Atom type parameter not found for '%s', using default radius %g\n", *atype[0], ATOM_DEFAULT_RADIUS);
      OUT.push_back( ATOM_DEFAULT_RADIUS );
    }
  }
  return OUT;
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

  // Get first ion
  DataSet_Coords* Ion1 = GetIonUnit( ion1name_, DSL );
  if (Ion1 == 0) {
    mprinterr("Error: Ion with name '%s' not found.\n", ion1name_.c_str());
    return 1;
  }
  // Check that there is a net charge
  double chargeIon1 = Ion1->Top().TotalCharge();
  if ( fabs(chargeIon1) < Constants::SMALL) {
    mprinterr("Error: Ion unit 1 '%s' does not have a net charge\n", Ion1->legend());
    return 1;
  }

  // Get second ion
  DataSet_Coords* Ion2 = 0;
  double chargeIon2 = 0.0;
  if (!ion2name_.empty()) {
    Ion2 = GetIonUnit( ion2name_, DSL );
    if (Ion2 == 0) {
      mprinterr("Error: Ion with name '%s' not found.\n", ion2name_.c_str());
      return 1;
    }
    // Check that there is a net charge
    chargeIon2 = Ion2->Top().TotalCharge();
    if ( fabs(chargeIon2) < Constants::SMALL) {
      mprinterr("Error: Ion unit 2 '%s' does not have a net charge\n", Ion2->legend());
      return 1;
    }
  }

  // Check that ion can actually neutralize
  int iIon1 = Nion1_;
  if (Nion1_ < 1) {
    if ( (chargeIon1 < 0 && totalCharge < 0) ||
         (chargeIon1 > 0 && totalCharge > 0) )
    {
      mprinterr("Error: First ion and system charges have same sign (%g and %g); can't neutralize.\n",
                chargeIon1, totalCharge);
      return 1;
    }
    // Get the nearest integer number of ions needed to neutralize the system.
    iIon1 = (int)lrint( fabs(totalCharge) / fabs(chargeIon1) );
    mprintf("\tNumber of %s ions required to neutralize: %i\n", ion1name_.c_str(), iIon1);
  }

  // Get atom radius for each ion
  typedef std::vector<double> Darray;

  Darray ion1radii, ion2radii;

  ion1radii = GetAtomRadii( Ion1, set0.AT() );

  if (Ion2 != 0)
    ion2radii = GetAtomRadii( Ion2, set0.AT() );

  mprintf("DEBUG: Nsolvent = %i\n", topOut.Nsolvent());

  return 0;
}
