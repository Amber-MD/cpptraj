#include <cmath>
#include <stack>
#include "Action_LipidOrder.h"
#include "CpptrajStdio.h"

const unsigned int Action_LipidOrder::MAX_H_ = 3;

/// CONSTRUCTOR
Action_LipidOrder::Action_LipidOrder() : axis_(DX) {}

// Action_LipidOrder::Help()
void Action_LipidOrder::Help() const {
  mprintf("\t[<name>] [<mask>] [{x|y|z}] [out <file>]\n");
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
  outfile_ = init.DFL().AddDataFile( actionArgs.GetStringKey("out") );
  mask_.SetMaskString( actionArgs.GetMaskNext() );
  dsname_ = actionArgs.GetStringNext();
  masterDSL_ = init.DslPtr();

  mprintf("    LIPIDORDER:\n");
  mprintf("\tCalculating lipid order parameters for lipids in mask '%s'\n", mask_.MaskString());
  const char AXISSTRING[3] = { 'X', 'Y', 'Z' };
  mprintf("\tCalculating with respect to the %c axis.\n", AXISSTRING[axis_]);
  return Action::OK;
}

// Action_LipidOrder::FindChain()
int Action_LipidOrder::FindChain( Npair const& chain ) {
  for (unsigned int idx = 0; idx != Types_.size(); idx++)
    if ( Types_[idx] == chain ) {
      mprintf("DEBUG: Existing chain: %s %s\n", *(Types_[idx].first), *(Types_[idx].second));
      return idx;
    }
  // New chain type
  mprintf("DEBIG: New chain: %s %s\n", *(chain.first), *(chain.second));
  Types_.push_back( chain );
  Chains_.push_back( ChainType() );
  return Types_.size()-1;
}

// Action_LipidOrder::Setup()
/** Within atom selection, attempt to identify C-HX in lipids "south"
  * of carbonyl.
  */
Action::RetType Action_LipidOrder::Setup(ActionSetup& setup)
{
  if (setup.Top().SetupCharMask( mask_ )) return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprintf("Warning: No atoms selected.\n");
    return Action::SKIP;
  }
  // Clear existing sites, but not data.
  Sites_.clear();
  // Loop over all molecules.
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
              // Determine if this is a new or existing chain type
              int chainIdx = FindChain( Npair(setup.Top().Res(atom.ResNum()).Name(), atom.Name()) );
              ChainType& Chain = Chains_[chainIdx];
              // Starting at the bonded carbon follow the chain down
              mprintf("DEBUG: Lipid chain type %i (%zu atoms) starting at (but not including) %s\n",
                      chainIdx, Chains_[chainIdx].size(), setup.Top().TruncResAtomName(at).c_str());
              int position = 0;
              std::stack<int> nextAtom;
              nextAtom.push( C_idx );
              while (!nextAtom.empty()) {
                int current_at = nextAtom.top();
                nextAtom.pop();
                Visited[current_at - offset] = true;
                Atom const& curr_atm = setup.Top()[current_at];
                // Add carbon if it does not yet exist in chain.
                if (position >= (int)Chain.size())
                  Chain.resize( position+1 );
                CarbonData& Cdata = Chain[position];
                if (!Cdata.Init())
                  Cdata.SetName( curr_atm.Name() );
                else if (Cdata.Name() != curr_atm.Name())
                  mprintf("Warning: Atom name %s at position %i (#%i) does not match %s\n",
                          *(curr_atm.Name()), position+1, current_at+1, Cdata.name());
                // Add site
                Sites_.push_back( CarbonSite(current_at, chainIdx, position) );
                CarbonSite& site = Sites_.back();
                mprintf("\t\t%i %s %s %i %i\n", position,
                        setup.Top().TruncResAtomNameNum(current_at).c_str(),
                        Cdata.name(), current_at, chainIdx);
                // Loop over atoms bonded to this carbon, add hydrogens to site
                for (Atom::bond_iterator bnd = setup.Top()[current_at].bondbegin();
                                         bnd != setup.Top()[current_at].bondend(); ++bnd)
                {
                  if (setup.Top()[*bnd].Element() == Atom::HYDROGEN)
                    site.AddHindex( *bnd );
                  else if (!Visited[*bnd - offset] &&
                           setup.Top()[*bnd].Element() == Atom::CARBON)
                    nextAtom.push( *bnd );
                }
                // Update number of hydrogens in carbon data
                Cdata.SetNumH( site.NumH() );
                position++;
              } // END loop over atom number stack
            } // END carbon is carbonyl
          } // END if carbon
        } // END if selected
      } // END loop over molecule atoms
    } // END mol not solvent and some of mol selected
  } // END loop over molecules

  mprintf("\t%zu chain types:\n", Types_.size());
  for (unsigned int idx = 0; idx != Types_.size(); idx++)
  {
    mprintf("\t[%u] %s %s (%zu)\n", idx, *(Types_[idx].first), *(Types_[idx].second),
            Chains_[idx].size());
    for (ChainType::const_iterator it = Chains_[idx].begin(); it != Chains_[idx].end(); ++it)
      mprintf("\t  %s\n", it->name());
  }
  return Action::OK;
}

// Action_LipidOrder::DoAction()
Action::RetType Action_LipidOrder::DoAction(int frameNum, ActionFrame& frm)
{
  // Loop over all carbon sites
  for (Carray::const_iterator site = Sites_.begin(); site != Sites_.end(); ++site)
  {
    Vec3 Cvec( frm.Frm().XYZ( site->Cidx() ) );
    CarbonData& cdata = Chains_[site->ChainIdx()][site->Position()];
    // Loop over hydrogens
    for (unsigned int i = 0; i != site->NumH(); i++)
    {
      // C-H unit vector
      Vec3 sx = Vec3(frm.Frm().XYZ( site->Hidx(i) )) - Cvec;
      sx.Normalize();
//      mprintf("DBG: %8i %8i %8.3f\n",site->Cidx(), site->Hidx(i), sx[axis_]);
      cdata.UpdateAngle(i, 0.5 * (3.0 * sx[axis_] * sx[axis_] - 1.0));
    }
    cdata.UpdateNvals();
  }
  return Action::OK;
}

//  Action_LipidOrder::Print()
void Action_LipidOrder::Print() {
  mprintf("    LIPIDORDER:\n");
  for (unsigned int idx = 0; idx != Types_.size(); idx++)
  {
    const char* resName = *(Types_[idx].first);
    const char* atmName = *(Types_[idx].second);
    // Create data set for chain type. Each hydrogen gets its own set.
    if (dsname_.empty())
      dsname_ = masterDSL_->GenerateDefaultName("LIPID");
    DataSet* DS[3];
    DS[0] = masterDSL_->AddSet(DataSet::DOUBLE, MetaData(dsname_, "H1", idx));
    DS[1] = masterDSL_->AddSet(DataSet::DOUBLE, MetaData(dsname_, "H2", idx));
    DS[2] = masterDSL_->AddSet(DataSet::DOUBLE, MetaData(dsname_, "H3", idx));
    if (DS[0] == 0 || DS[1] == 0 || DS[2] == 0) {
      mprinterr("Error: Could not create data sets for chain %s %s\n", resName, atmName);
      continue;
    }
    if (outfile_ != 0) {
      outfile_->AddDataSet( DS[0] );
      outfile_->AddDataSet( DS[1] );
      outfile_->AddDataSet( DS[2] );
    }
    std::string prefix = Types_[idx].first.Truncated() + "_" + Types_[idx].second.Truncated();
    DS[0]->SetLegend( prefix + "_H1" );
    DS[1]->SetLegend( prefix + "_H2" );
    DS[2]->SetLegend( prefix + "_H3" );
    // Loop over carbons in chain
    int pos = 0;
    for (ChainType::const_iterator it = Chains_[idx].begin(); it != Chains_[idx].end(); ++it)
    {
      mprintf("\t%s %s (%zu)", resName, it->name(), it->Nvals());
      if (it->Nvals() > 0) {
        double avg, stdev;
        for (unsigned int i = 0; i != MAX_H_; i++) {
          //mprintf(" %15s", setup.Top().TruncResAtomName(site->Hidx(i)).c_str());
          if ( i < it->NumH() )
            avg = it->Avg(i, stdev);
          mprintf("  %10.7f %10.7f  ", avg, stdev);
          DS[i]->Add(pos, &avg);
        }
        mprintf("\n");
      }
    }
  }
}


// =============================================================================
/// CONSTRUCTOR
Action_LipidOrder::CarbonSite::CarbonSite() : c_idx_(-1), chainIdx_(-1), position_(-1), nH_(0)
{
  std::fill(h_idx_, h_idx_ + MAX_H_, -1);
}

/** CONSTRUCTOR - take carbon index */
Action_LipidOrder::CarbonSite::CarbonSite(int c, int i, int p) :
  c_idx_(c), chainIdx_(i), position_(p), nH_(0)
{
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
/// CONSTRUCTOR
Action_LipidOrder::CarbonData::CarbonData() :
  nvals_(0), nH_(0), init_(false)
{
  std::fill( sum_,  sum_  + MAX_H_, 0.0 );
  std::fill( sum2_, sum2_ + MAX_H_, 0.0 );
}

/** \return average and calc standard dev. for specified hydrogen */
double Action_LipidOrder::CarbonData::Avg(int i, double& stdev) const {
  double avg = sum_[i] / (double)nvals_;
  stdev = sum2_[i] / (double)nvals_;
  stdev -= (avg * avg);
  if (stdev > 0.0)
    stdev = sqrt( stdev );
  else
    stdev = 0.0;
  return avg;
}
