#include "Action_LipidOrder.h"
#include "CpptrajStdio.h"

// Action_LipidOrder::Help()
void Action_LipidOrder::Help() const {

}

// Action_LipidOrder::Init()
Action::RetType Action_LipidOrder::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  mask_.SetMaskString( actionArgs.GetMaskNext() );

  mprintf("    LIPIDORDER:\n");
  return Action::OK;
}

// Action_LipidOrder::Setup()
Action::RetType Action_LipidOrder::Setup(ActionSetup& setup)
{
  // Within selection, attempt to identify C-H in lipids
  // "south" of the carbonyl
  if (setup.Top().SetupCharMask( mask_ )) return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprintf("Warning: No atoms selected.\n");
    return Action::SKIP;
  }
  for (Topology::mol_iterator mol = setup.Top().MolStart();
                              mol != setup.Top().MolEnd(); ++mol)
  {
    if (!mol->IsSolvent() && mask_.AtomsInCharMask(mol->BeginAtom(), mol->EndAtom()))
    {
      std::vector<bool> Visited(mol->NumAtoms(), false);
      int idx = 0;
      for (int at = mol->BeginAtom(); at != mol->EndAtom(); ++at, ++idx)
      {
        Visited[idx] = true;
        if (mask_.AtomInCharMask( at )) {
          Atom const& atom = setup.Top()[at];
          if (atom.Element() == Atom::CARBON) {
            // Should be bonded to two oxygens and 1 carbon
            int n_O = 0;
            int n_C = 0;
            int C_idx = 0;
            for (Atom::bond_iterator bnd = atom.bondbegin();
                                     bnd != atom.bondend(); ++bnd)
            {
              if (setup.Top()[*bnd].Element() == Atom::OXYGEN)
                n_O++;
              else if (setup.Top()[*bnd].Element() == Atom::CARBON) {
                n_C++;
                C_idx = *bnd;
              } else {
                n_O = 0;
                n_C = 0;
                break;
              }
            }
            if (n_O == 2 && n_C == 1) {
              // Starting at the carbon follow the chain down
              mprintf("DEBUG: Lipid chain starting at %s\n",
                      setup.Top().TruncResAtomName(C_idx).c_str());
            }
          } // END if Carbon
        } // END if selected
      } // END loop over molecule atoms
    } // END mol not solvent and some of mol selected
  } // END loop over molecules
  return Action::ERR;
}

// Action_LipidOrder::DoAction()
Action::RetType Action_LipidOrder::DoAction(int frameNum, ActionFrame& frm)
{
  return Action::ERR;
}
