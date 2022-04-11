#include "Exec_HmassRepartition.h"
#include "CharMask.h"
#include "CpptrajStdio.h"

// Exec_HmassRepartition::Help()
void Exec_HmassRepartition::Help() const
{
  mprintf("\t[%s]\n"
          "\t[<mask>] [hmass <hydrogen new mass>] [dowater]\n", DataSetList::TopIdxArgs);
  mprintf("  Perform hydrogen mass repartitioning for atoms selected by <mask>\n"
          "  (all solute atoms by default).\n");
}

// Exec_HmassRepartition::Execute()
Exec::RetType Exec_HmassRepartition::Execute(CpptrajState& State, ArgList& argIn)
{
  double hmass = argIn.getKeyDouble("hmass", 3.024);
  bool do_water = argIn.hasKey("dowater");

  std::string maskArg = argIn.GetMaskNext();
  CharMask mask;
  if (mask.SetMaskString(maskArg)) {
    mprinterr("Error: Could not process mask string.\n");
    return CpptrajState::ERR;
  }
  // Get Topology
  Topology* topIn = State.DSL().GetTopByIndex( argIn );
  if (topIn == 0) return CpptrajState::ERR;

  mprintf("\tPerforming hydrogen mask repartitioning for atoms selected by mask '%s'\n",
          mask.MaskString());
  mprintf("\tTopology is '%s'\n", topIn->c_str());
  mprintf("\tHydrogen masses will be changed to %f amu\n", hmass);
  if (do_water)
    mprintf("\tRepartitioning solvent hydrogen masses as well.\n");
  else
    mprintf("\tSkipping solvent.\n");

  if (topIn->SetupCharMask( mask )) {
    mprinterr("Error: Could not set up atom mask.\n");
    return CpptrajState::ERR;
  }
  mask.MaskInfo();
  if (mask.None()) {
    mprintf("Warning: No atoms selected.\n");
    return CpptrajState::OK;
  }

  if (repartition(*topIn, hmass, mask, do_water, State.Debug()))
    return CpptrajState::ERR; 

  return CpptrajState::OK;
}

/** Do hydrogen mass repartitioning. Change mass of all hydrogens to 
  * the given mass. Adjust mass of bonded heavy atom to maintain the
  * same overall mass.
  */
int Exec_HmassRepartition::repartition(Topology& topIn, double hmass, CharMask const& cmaskIn,
                                       bool do_water, int debugIn)
{
  AtomMask amask( cmaskIn.ConvertToIntMask(), cmaskIn.Natom() );

  unsigned int n_h_changed = 0;
  unsigned int n_heavy_changed = 0;
  for (AtomMask::const_iterator at = amask.begin(); at != amask.end(); ++at)
  {
    // Skip solvent
    if (!do_water) {
      int molnum = topIn[*at].MolNum();
      if (topIn.Mol(molnum).IsSolvent()) continue;
    }

    if (topIn[*at].Element() != Atom::HYDROGEN) {
      Atom& heavyAtom = topIn.SetAtom(*at);
      double delta = 0;
      int n_selected_h_atoms = 0;
      for (Atom::bond_iterator bat = heavyAtom.bondbegin();
                               bat != heavyAtom.bondend(); ++bat)
      {
        if (cmaskIn.AtomInCharMask(*bat) && topIn[*bat].Element() == Atom::HYDROGEN)
        {
          n_selected_h_atoms++;
          Atom& hAtom = topIn.SetAtom(*bat);
          double diff = hmass - hAtom.Mass();
          hAtom.SetMass( hmass );
          n_h_changed++;
          delta += diff;
        }
      }
      if (n_selected_h_atoms > 0) {
        if (debugIn > 0)
          mprintf("\tHeavy atom %s bonded to %i hydrogens, subtract %f amu.\n",
                  topIn.TruncResAtomNameNum(*at).c_str(),
                  n_selected_h_atoms, delta);
        double newMass = heavyAtom.Mass() - delta;
        // Sanity check
        if (newMass < 0) {
          mprinterr("Error: New mass for '%s' would be less than 0.\n",
                    topIn.TruncResAtomNameNum(*at).c_str());
          return 1;
        }
        heavyAtom.SetMass( newMass );
        n_heavy_changed++;
      }
    }
  }
  mprintf("\t%u hydrogen masses modified.\n", n_h_changed);
  mprintf("\t%u heavy atom masses modified.\n", n_heavy_changed);
  
  return 0;
}
