#include "AddIons.h"
#include "../Parm/ParameterSet.h"
#include "../Parm/ParmHolder.h"
#include "../AtomType.h"
#include "../Constants.h"
#include "../CpptrajStdio.h"
#include "../DataSet_Coords.h"
#include "../DataSetList.h"
#include "../DistRoutines.h"
#include "../Topology.h"
#include <cmath> //fabs, lrint
#include <cstdlib> // FIXME debug

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
AddIons::AddIons() :
  Nion1_(0),
  Nion2_(0),
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
    Nion2_ = 0;
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
  if (separation_ > 0.0)
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
  // Check that there are solvent molecules.
  int nsolvent = topOut.Nsolvent();
  if (nsolvent < 1) {
    mprinterr("Error: No solvent present in '%s'. Add solvent first.\n", topOut.c_str());
    return 1;
  }
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
  int iIon2 = Nion2_; 
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

  // Check that there is enough solvent to swap
  int nIons = iIon1 + iIon2;
  if ( nIons == nsolvent )
    mprintf("Warning: # of ions to add (%i) is same as number of solvent molecules (%i)\n", nIons, nsolvent);
  else if (nIons > nsolvent) {
    mprinterr("Error: # of ions to add (%i) is larger than number of solvent molecules (%i)\n", nIons, nsolvent);
    return 1;
  }

  // Get atom radius for each ion
  typedef std::vector<double> Darray;

  Darray ion1radii, ion2radii;

  ion1radii = GetAtomRadii( Ion1, set0.AT() );

  if (Ion2 != 0)
    ion2radii = GetAtomRadii( Ion2, set0.AT() );

  mprintf("DEBUG: Nsolvent = %i\n", topOut.Nsolvent());

  std::vector<int> solventMolNums = topOut.SolventMolNums();
  // DEBUG
  mprintf("DEBUG: Solvent molecule #s:");
  for (std::vector<int>::const_iterator it = solventMolNums.begin(); it != solventMolNums.end(); ++it)
    mprintf(" %i", *it);
  mprintf("\n"); // DEBUG

  mprintf("Adding %d counter ions to \"%s\". %d solvent molecules will remain.\n",
          iIon1 + iIon2, topOut.c_str(), (int)topOut.SolventMolNums().size() - iIon1 - iIon2);

  typedef std::vector<Vec3> Varray;
  Varray ionPositions;
  if (separation_ > 0.0)
    ionPositions.reserve( nIons );
  double cut2 = separation_ * separation_;

  // FIXME DEBUG
  srand( 1 );

  int failCounter = 0;
  while ( iIon1 || iIon2 ) {
    if (iIon1) {
      if (place_ion(iIon1, failCounter, ionPositions, Ion1, topOut, frameOut, solventMolNums, cut2, nIons))
        return 1;
    }
    if (iIon2) {
      if (place_ion(iIon2, failCounter, ionPositions, Ion2, topOut, frameOut, solventMolNums, cut2, nIons))
        return 1;
    }
  }

  return 0;
}

/** Select a solvent molecule to swap with that is not too close to other ions positions. */
int AddIons::place_ion(int& iIon1, int& failCounter, Varray& ionPositions,
                       DataSet_Coords* ionCrd, Topology& topOut, Frame& frameOut,
                       std::vector<int>& solventMolNums, double cut2, int nIons)
const
{
  // Pick a random solvent molecule to replace
  int solvIdx = RNG_.rn_num() % (int)solventMolNums.size();
  int ntries = 0;
  while ( solventMolNums[solvIdx] == -1 ) {
    ntries++;
    if (ntries > 100) {
      mprinterr("Error: Could not find a solvent molecule to swap with after 100 tries.\n");
      return 1;
    }
    solvIdx = RNG_.rn_num() % (int)solventMolNums.size();
  }
  
  // Do not try it again. It will either be used or be invalid due to distances.
  int solventMolNum = solventMolNums[solvIdx];
  mprintf("DEBUG: ntries=%i Trying swap %s with solvent molecule %i\n", ntries, ionCrd->legend(), solventMolNum);
  solventMolNums[solvIdx] = -1;
 
  // Get position of solvent residue atom. TODO should this be the geometric center?
  Molecule const& solventMol = topOut.Mol( solventMolNum );
  int firstAt = solventMol.MolUnit().Front();
  Vec3 pos( frameOut.XYZ(firstAt) );

  // Check that this position isnt too close to other positions
  bool placeIon = true;
  for (Varray::const_iterator previousPos = ionPositions.begin(); previousPos != ionPositions.end(); ++previousPos)
  {
    // TODO imaging
    double dist2 = DIST2_NoImage( pos, *previousPos );
    if (dist2 < cut2) {
      ++failCounter;
      placeIon = false;
      break;
    }
  }

  if (placeIon) {
    std::vector<int> validIonResidues;
    validIonResidues.reserve( ionCrd->Top().Nres() );
    for (int ires = 0; ires != ionCrd->Top().Nres() ; ++ires)
      validIonResidues.push_back( ires );
    Frame ionFrame = ionCrd->AllocateFrame(); //TODO allocate outside this routine
    ionCrd->GetFrame(0, ionFrame);
    ionFrame.CenterOnPoint( pos, false ); // Use geometric center
    Vec3 debugVec = ionFrame.VGeometricCenter();
    debugVec.Print("DEBUG: check ion center");
    mprintf("%zu: Placed %s in %s at (%4.2lf, %4.2lf, %4.2lf).\n", iIon1, ionCrd->legend(), topOut.c_str(),
            pos[0], pos[1], pos[2]);
    if (separation_ > 0.0)
      ionPositions.push_back( pos );
    topOut.AddResidues( ionCrd->Top(), validIonResidues, frameOut, ionFrame, false );
    --iIon1;
  }

  // Safety valve
  if (failCounter > 100) {
    mprinterr("Error: Could not place %i ions with a minimum separation of %f after 100 tries.\n", nIons, separation_);
    return 1;
  }
  
  return 0;
}
