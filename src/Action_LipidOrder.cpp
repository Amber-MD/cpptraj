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

//  Action_LipidOrder::FollowChain()
void Action_LipidOrder::FollowChain(int idx, Topology const& top, std::vector<bool> Visited,
                                    int offset)
{
  Visited[idx - offset] = true;
  Atom const& atom = top[idx];
  // Get site
  Cmap::iterator it = Carbons_.lower_bound( atom.Name() );
  if (it == Carbons_.end() || it->first != atom.Name())
    // New site 
    it = Carbons_.insert( it, Cpair(atom.Name(), CarbonList(top.Res(atom.ResNum()).Name())) );

  // Add this site
  mprintf("\t\t%s\n", top.TruncResAtomName(idx).c_str());
  CarbonSite& site = it->second.AddSite( idx );
  // Loop over bonded atoms 
  for (Atom::bond_iterator bnd = top[idx].bondbegin();
                           bnd != top[idx].bondend(); ++bnd)
  {
    if (top[*bnd].Element() == Atom::HYDROGEN)
      site.AddHindex( *bnd );
    else if (!Visited[*bnd - offset] && top[*bnd].Element() == Atom::CARBON)
      FollowChain(*bnd, top, Visited, offset);
  }
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
      mprintf("  Mol %zu\n", mol - setup.Top().MolStart() + 1);
      std::vector<bool> Visited(mol->NumAtoms(), false);
      int offset = mol->BeginAtom();
      // First mark any atoms not in the mask as visited.
      for (int at = mol->BeginAtom(); at != mol->EndAtom(); ++at)
        if (!mask_.AtomInCharMask( at ))
          Visited[at - offset] = true;
      // Search for the start of lipid chains
      for (int at = mol->BeginAtom(); at != mol->EndAtom(); ++at)
      {
        if (!Visited[at - offset]) {
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
              Visited[at - offset] = true;
              // Starting at the bonded carbon follow the chain down
              mprintf("DEBUG: Lipid chain starting at %s\n",
                      setup.Top().TruncResAtomName(C_idx).c_str());
              FollowChain(C_idx, setup.Top(), Visited, offset);
            }
          } // END if Carbon
        } // END if selected
      } // END loop over molecule atoms
    } // END mol not solvent and some of mol selected
  } // END loop over molecules

  mprintf("\t%zu carbon types:\n", Carbons_.size());
  for (Cmap::const_iterator it = Carbons_.begin(); it != Carbons_.end(); ++it)
  {
    mprintf("\t%s %s (%zu)\n", it->second.resName(), *(it->first), it->second.Nsites());
    for (Carray::const_iterator site = it->second.begin(); site != it->second.end(); ++site)
    {
      mprintf("\t%15s", setup.Top().TruncResAtomName(site->Cidx()).c_str());
      if (site->H1idx() != -1) mprintf(" %15s", setup.Top().TruncResAtomName(site->H1idx()).c_str());
      if (site->H2idx() != -1) mprintf(" %15s", setup.Top().TruncResAtomName(site->H2idx()).c_str());
      if (site->H3idx() != -1) mprintf(" %15s", setup.Top().TruncResAtomName(site->H3idx()).c_str());
      mprintf("\n");
    }
  }
  return Action::ERR;
}

// Action_LipidOrder::DoAction()
Action::RetType Action_LipidOrder::DoAction(int frameNum, ActionFrame& frm)
{
  return Action::ERR;
}
