#include <cmath>
#include "Action_LipidOrder.h"
#include "CpptrajStdio.h"

const unsigned int Action_LipidOrder::MAX_H_ = 3;

Action_LipidOrder::Action_LipidOrder() : axis_(DX) {}

// Action_LipidOrder::Help()
void Action_LipidOrder::Help() const {

}

// Action_LipidOrder::Init()
Action::RetType Action_LipidOrder::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  if (actionArgs.hasKey("x") )
    axis_ = DX;
  else if (actionArgs.hasKey("y") )
    axis_ = DY;
  else if (actionArgs.hasKey("z") )
    axis_ = DZ;
  else
    axis_ = DZ;
  mask_.SetMaskString( actionArgs.GetMaskNext() );

  mprintf("    LIPIDORDER:\n");
  mprintf("\tCalculating lipid order parameters for lipids in mask '%s'\n", mask_.MaskString());
  const char AXISSTRING[3] = { 'X', 'Y', 'Z' };
  mprintf("\tCalculating with respect to the %c axis.\n", AXISSTRING[axis_]);
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
  mprintf("\t\t%s\n", top.TruncResAtomNameNum(idx).c_str());
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
            if (n_O == 2 && n_C == 1 && mask_.AtomInCharMask(C_idx)) {
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
      mprintf("\t%20s", setup.Top().TruncResAtomNameNum(site->Cidx()).c_str());
      for (unsigned int i = 0; i != site->NumH(); i++)
        mprintf(" %20s", setup.Top().TruncResAtomNameNum(site->Hidx(i)).c_str());
      mprintf("\n");
    }
  }
  return Action::OK;
}

// Action_LipidOrder::DoAction()
Action::RetType Action_LipidOrder::DoAction(int frameNum, ActionFrame& frm)
{
  // Loop over types
  for (Cmap::iterator it = Carbons_.begin(); it != Carbons_.end(); ++it)
  {
    // Loop over sites
    for (Carray::iterator site = it->second.begin(); site != it->second.end(); ++site)
    {
      Vec3 Cvec( frm.Frm().XYZ( site->Cidx() ) );
      // Loop over hydrogens
      for (unsigned int i = 0; i != site->NumH(); i++)
      {
        // C-H unit vector
        Vec3 sx = Cvec - Vec3(frm.Frm().XYZ( site->Hidx(i) ));
        sx.Normalize();
        mprintf("DBG: %8i %8i %8.3f\n",site->Cidx(), site->Hidx(i), sx[axis_]); 
        it->second.UpdateAngle(i, 0.5 * (3.0 * sx[axis_] * sx[axis_] - 1.0));
      }
      it->second.UpdateNvals();
    }
  }
  return Action::ERR;
}

//  Action_LipidOrder::Print()
void Action_LipidOrder::Print() {
  mprintf("    LIPIDORDER:\n");
  for (Cmap::const_iterator it = Carbons_.begin(); it != Carbons_.end(); ++it)
  {
    mprintf("\t%s %s (%zu)", it->second.resName(), *(it->first), it->second.Nvals());
    if (it->second.Nvals() > 0) {
      Carray::const_iterator first = it->second.begin();
      // FIXME assuming all sites have same # H
      for (unsigned int i = 0; i != first->NumH(); i++) {
        //mprintf(" %15s", setup.Top().TruncResAtomName(site->Hidx(i)).c_str());
        double avg, stdev;
        avg = it->second.Avg(i, stdev);
        mprintf("  %10.7f %10.7f  ", avg, stdev);
      }
      mprintf("\n");
    }
  }
}


// =============================================================================

/// CONSTRUCTOR
Action_LipidOrder::CarbonSite::CarbonSite() : c_idx_(-1), nH_(0) {
  std::fill(h_idx_, h_idx_ + MAX_H_, -1);
}

/** CONSTRUCTOR - take carbon index */
Action_LipidOrder::CarbonSite::CarbonSite(int c) : c_idx_(c), nH_(0) {
  std::fill(h_idx_, h_idx_ + MAX_H_, -1);
}

/** Add hydrogen index */
void Action_LipidOrder::CarbonSite::AddHindex(int h) {
  h_idx_[nH_] = h;
  if (nH_ < 3) nH_++;
  else {
    mprintf("Warning: Attempting to add 4th hydrogen (index %i) to carbon index %i\n",
            h+1, c_idx_+1);
    mprintf("Warning: Replaced existing index.\n");
  }
}

// =============================================================================
// CONSTRUCTOR
Action_LipidOrder::CarbonList::CarbonList() : nvals_(0) {
  std::fill(sum_, sum_ + MAX_H_, 0.0);
  std::fill(sum2_, sum2_ + MAX_H_, 0.0);
}

/** CONSTRUCTOR - Take residue name associated with carbon site. */
Action_LipidOrder::CarbonList::CarbonList(NameType const& n) : resname_(n), nvals_(0)
{
  std::fill(sum_, sum_ + MAX_H_, 0.0);
  std::fill(sum2_, sum2_ + MAX_H_, 0.0);
}

/** \return average for hydrogen */
double Action_LipidOrder::CarbonList::Avg(int i, double& stdev) const {
  double avg = sum_[i] / (double)nvals_;
  stdev = sum2_[i] / (double)nvals_;
  stdev -= (avg * avg);
  if (stdev > 0.0)
    stdev = sqrt( stdev );
  else
    stdev = 0.0;
  return avg;
}

// Action_LipidOrder::CarbonList::AddSite()
Action_LipidOrder::CarbonSite& Action_LipidOrder::CarbonList::AddSite(int cidx) {
  // Does this site already exist?
  for (Carray::const_iterator site = sites_.begin(); site != sites_.end(); ++site)
    if (site->Cidx() == cidx) return sites_[site - sites_.begin()];
  // Does not exist, add.
  sites_.push_back( CarbonSite(cidx) );
  return sites_.back();
}
