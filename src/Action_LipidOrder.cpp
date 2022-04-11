#include <cmath>
#include <stack>
#include "Action_LipidOrder.h"
#include "CpptrajStdio.h"
#ifdef _OPENMP
# include <omp.h>
#endif

const unsigned int Action_LipidOrder::MAX_H_ = 3;

/// CONSTRUCTOR
Action_LipidOrder::Action_LipidOrder() :
  axis_(DX),
  outfile_(0),
  debug_(0),
  report_p2_(false),
  nthreads_(1)
{}

// Action_LipidOrder::Help()
void Action_LipidOrder::Help() const {
  mprintf("\t[<name>] [<mask>] [{x|y|z}] [out <file>] [p2]\n"
          "  Calculate lipid order parameters SCD (|<P2>|) for lipid chains in mask\n"
          "  <mask>. Lipid chains are identified by carboxyl groups, i.e.\n"
          "  O-(C=O)-C1-...-CN, where C1 is the first carbon in the acyl chain\n"
          "  and CN is the last. Order parameters will be determined for each\n"
          "  hydrogen bonded to each carbon. If 'p2' is specified the raw <P2>\n"
          "  values will be reported.\n");
}

// Action_LipidOrder::Init()
Action::RetType Action_LipidOrder::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  masterDSL_ = init.DslPtr();
# ifdef MPI
  trajComm_ = init.TrajComm();
# endif
  if (actionArgs.hasKey("x") )
    axis_ = DX;
  else if (actionArgs.hasKey("y") )
    axis_ = DY;
  else if (actionArgs.hasKey("z") )
    axis_ = DZ;
  else
    axis_ = DZ;
  report_p2_ = actionArgs.hasKey("p2");
  outfile_ = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  if (mask_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;
  dsname_ = actionArgs.GetStringNext();

# ifdef _OPENMP
  // Each thread needs space to calc scd each frame 
# pragma omp parallel
  {
# pragma omp master
  {
  nthreads_ = omp_get_num_threads();
  }
  }
# endif

  mprintf("    LIPIDORDER:\n");
  const char* typestr;
  if (report_p2_)
    typestr = "<P2>";
  else
    typestr = "SCD=|<P2>|";
  mprintf("\tCalculating lipid order parameters (%s) for lipids in mask '%s'\n",
          typestr, mask_.MaskString());
  const char AXISSTRING[3] = { 'X', 'Y', 'Z' };
  mprintf("\tCalculating with respect to the %c axis.\n", AXISSTRING[axis_]);
  if (!dsname_.empty())
    mprintf("\tData saved in sets named '%s'\n", dsname_.c_str());
  if (outfile_ != 0)
    mprintf("\tOutput to file '%s'\n", outfile_->DataFilename().full());
  if (nthreads_ > 1)
    mprintf("\tParallelizing calculation with %i threads.\n", nthreads_);

  return Action::OK;
}

// Action_LipidOrder::FindChain()
int Action_LipidOrder::FindChain( Npair const& chain ) {
  for (unsigned int idx = 0; idx != Types_.size(); idx++)
    if ( Types_[idx] == chain ) {
      if (debug_ > 0)
        mprintf("DEBUG: Existing chain: %s %s\n", *(Types_[idx].first), *(Types_[idx].second));
      Nchains_[idx]++;
      return idx;
    }
  // New chain type
  if (debug_ > 0)
    mprintf("DEBUG: New chain: %s %s\n", *(chain.first), *(chain.second));
  Types_.push_back( chain );
  Chains_.push_back( ChainType() );
  Nchains_.push_back( 1 );
  return Types_.size()-1;
}

// Action_LipidOrder::Setup()
/** Within atom selection, attempt to identify C-HX in lipids "south"
  * of carboxyl.
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
  Nchains_.clear();
  // Loop over all molecules.
  unsigned int nChains = 0;
  for (Topology::mol_iterator mol = setup.Top().MolStart();
                              mol != setup.Top().MolEnd(); ++mol)
  {
    if (!mol->IsSolvent() && mask_.AtomsInCharMask(mol->MolUnit()))
    {
      if (debug_ > 0)
        mprintf("  Mol %zu\n", mol - setup.Top().MolStart() + 1);
      // Array to keep track of visited atoms. Since molecules do not need
      // to be contiguous have the array correspond to all atoms. Initially
      // set it up so that all atoms are visited.
      std::vector<bool> Visited(setup.Top().Natom(), true);
      // Mark any atoms in the mask as not visited.
      for (Unit::const_iterator seg = mol->MolUnit().segBegin();
                                seg != mol->MolUnit().segEnd(); ++seg)
      {
        for (int at = seg->Begin(); at != seg->End(); ++at)
          if (mask_.AtomInCharMask( at ))
            Visited[at] = false;
      }
      // Search for the start of lipid chains
      for (Unit::const_iterator seg = mol->MolUnit().segBegin();
                                seg != mol->MolUnit().segEnd(); ++seg)
      {
        for (int at = seg->Begin(); at != seg->End(); ++at)
        {
          if (!Visited[at]) {
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
                nChains++;
                Visited[at] = true;
                // Determine if this is a new or existing chain type by first
                // carbon residue name and carboxyl carbon atom name.
                int cresnum = setup.Top()[C_idx].ResNum();
                int chainIdx = FindChain( Npair(setup.Top().Res(cresnum).Name(), atom.Name()) );
                ChainType& Chain = Chains_[chainIdx];
                // Starting at the bonded carbon follow the chain down
                if (debug_ > 0)
                  mprintf("DEBUG: Lipid chain type %i (%zu atoms) starting at"
                          " (but not including) %s\n", chainIdx, Chains_[chainIdx].size(),
                          setup.Top().TruncResAtomName(at).c_str());
                int position = 0;
                std::stack<int> nextAtom;
                nextAtom.push( C_idx );
                while (!nextAtom.empty()) {
                  int current_at = nextAtom.top();
                  nextAtom.pop();
                  Visited[current_at] = true;
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
                  if (debug_ > 1)
                    mprintf("\t\t%i %s %s %i %i\n", position,
                            setup.Top().TruncResAtomNameNum(current_at).c_str(),
                            Cdata.name(), current_at, chainIdx);
                  // Loop over atoms bonded to this carbon, add hydrogens to site
                  for (Atom::bond_iterator bnd = setup.Top()[current_at].bondbegin();
                                           bnd != setup.Top()[current_at].bondend(); ++bnd)
                  {
                    if (setup.Top()[*bnd].Element() == Atom::HYDROGEN)
                      site.AddHindex( *bnd );
                    else if (!Visited[*bnd] &&
                             setup.Top()[*bnd].Element() == Atom::CARBON)
                      nextAtom.push( *bnd );
                  }
                  // Update number of hydrogens in carbon data
                  Cdata.SetNumH( site.NumH() );
                  position++;
                } // END loop over atom number stack
              } // END carbon is carboxyl
            } // END if carbon
          } // END if selected
        } // END loop over molecule segment atoms
      } // END loop over molecule segments
    } // END mol not solvent and some of mol selected
  } // END loop over molecules

  mprintf("\t%u chains.\n", nChains);
  mprintf("\t%zu chain types:\n", Types_.size());
  for (unsigned int idx = 0; idx != Types_.size(); idx++)
  {
    mprintf("\t[%u] %i chains, Residue %s, carboxyl atom %s (%zu carbons)\n", idx,
            Nchains_[idx],
            *(Types_[idx].first),
            *(Types_[idx].second), Chains_[idx].size());
    mprintf("\t  %-4s %-4s %2s\n", "Pos.", "Name", "#H");
    unsigned int pos = 2;
    for (ChainType::const_iterator it = Chains_[idx].begin(); it != Chains_[idx].end(); ++it, ++pos)
      mprintf("\t  %-4u %-4s %2u\n", pos, it->name(), it->NumH());
  }

  return Action::OK;
}

// Action_LipidOrder::DoAction()
Action::RetType Action_LipidOrder::DoAction(int frameNum, ActionFrame& frm)
{
  // Zero all temp arrays.
  for (ChainArray::iterator chn = Chains_.begin(); chn != Chains_.end(); ++chn)
    for (ChainType::iterator it = chn->begin(); it != chn->end(); ++it)
#     ifdef _OPENMP
      it->AllocateTempSpace(nthreads_);
#     else
      it->AllocateTempSpace();
#     endif

  // Loop over all carbon sites
  int idx;
  int maxIdx = (int)Sites_.size();
  int offset = 0;
# ifdef _OPENMP
  int mythread;
# pragma omp parallel private(mythread, offset, idx)
  {
    mythread = omp_get_thread_num();
#   pragma omp for
# endif
    for (idx = 0; idx < maxIdx; idx++)
    {
      CarbonSite const& site = Sites_[idx];
      Vec3 Cvec( frm.Frm().XYZ( site.Cidx() ) );
      CarbonData& cdata = Chains_[site.ChainIdx()][site.Position()];
#     ifdef _OPENMP
      offset = mythread * cdata.NumH();
#     endif
      // Loop over hydrogens
      for (unsigned int i = 0; i != site.NumH(); i++)
      {
        // C-H unit vector
        Vec3 sx = Vec3(frm.Frm().XYZ( site.Hidx(i) )) - Cvec; 
        sx.Normalize();
        //cdata.UpdateAngle(offset+i, 0.5 * (3.0 * sx[axis_] * sx[axis_] - 1.0));
        double scdval = 0.5 * (3.0 * sx[axis_] * sx[axis_] - 1.0);
        cdata.IncrementTempBy( offset+i, scdval );
      }
      //cdata.UpdateNvals( mythread );
    }
# ifdef _OPENMP
  } // END pragma omp parallel
# endif

  // Normalize by number of chains
  for (unsigned int cidx = 0; cidx != Chains_.size(); cidx++)
  {
    double dval = 1.0 / (double)Nchains_[cidx];
    for (unsigned int pos = 0; pos != Chains_[cidx].size(); pos++)
    {
      CarbonData& cdata = Chains_[cidx][pos];
#     ifdef _OPENMP
      cdata.ConsolidateTemp( nthreads_ );
#     endif
      for (unsigned int i = 0; i != cdata.NumH(); i++)
        cdata.UpdateAngle(i, cdata.TempSCD(i) * dval);
      cdata.UpdateNvals();
    }
  }

  return Action::OK;
}

#ifdef MPI
// Action_LipidOrder::SyncAction()
int Action_LipidOrder::SyncAction() {
  int total = 3;
  double dbuf[3];
  unsigned int ubuf[1];
  for (unsigned int idx = 0; idx != Types_.size(); idx++)
  {
    for (ChainType::iterator it = Chains_[idx].begin(); it != Chains_[idx].end(); ++it)
    {
      trajComm_.ReduceMaster( dbuf, it->Sptr(), total, MPI_DOUBLE, MPI_SUM );
      if (trajComm_.Master()) std::copy( dbuf, dbuf+total, it->Sptr() );
      trajComm_.ReduceMaster( dbuf, it->S2ptr(), total, MPI_DOUBLE, MPI_SUM );
      if (trajComm_.Master()) std::copy( dbuf, dbuf+total, it->S2ptr() );
      trajComm_.ReduceMaster( ubuf, it->Nptr(), 1, MPI_UNSIGNED, MPI_SUM );
      if (trajComm_.Master()) std::copy( ubuf, ubuf+1, it->Nptr() );
    }
  }
  return 0;
}
#endif /* MPI */

//  Action_LipidOrder::Print()
void Action_LipidOrder::Print() {
  if (debug_ > 0) mprintf("    LIPIDORDER:\n");
  Dimension Xaxis(2, 1, "Cn");
  for (unsigned int idx = 0; idx != Types_.size(); idx++)
  {
    const char* resName = *(Types_[idx].first);
    const char* atmName = *(Types_[idx].second);
    // Create data set for chain type. Each hydrogen gets its own set.
    DataSet* DS[3];
    DataSet* SD[3];
    if (dsname_.empty())
      dsname_ = masterDSL_->GenerateDefaultName("LIPID");
    DS[0] = masterDSL_->AddSet(DataSet::DOUBLE, MetaData(dsname_, "H1", idx));
    SD[0] = masterDSL_->AddSet(DataSet::DOUBLE, MetaData(dsname_, "SDH1", idx));
    DS[1] = masterDSL_->AddSet(DataSet::DOUBLE, MetaData(dsname_, "H2", idx));
    SD[1] = masterDSL_->AddSet(DataSet::DOUBLE, MetaData(dsname_, "SDH2", idx));
    DS[2] = masterDSL_->AddSet(DataSet::DOUBLE, MetaData(dsname_, "H3", idx));
    SD[2] = masterDSL_->AddSet(DataSet::DOUBLE, MetaData(dsname_, "SDH3", idx));
    if (DS[0] == 0 || DS[1] == 0 || DS[2] == 0 || SD[0] == 0 || SD[1] == 0 || SD[2] == 0) {
      mprinterr("Error: Could not create data sets for chain %s %s\n", resName, atmName);
      continue;
    }
    if (outfile_ != 0) {
      outfile_->AddDataSet( DS[0] );
      outfile_->AddDataSet( SD[0] );
      outfile_->AddDataSet( DS[1] );
      outfile_->AddDataSet( SD[1] );
      outfile_->AddDataSet( DS[2] );
      outfile_->AddDataSet( SD[2] );
    }
    std::string prefix = Types_[idx].first.Truncated() + "_" + Types_[idx].second.Truncated();
    DS[0]->SetLegend( prefix + "_H1" );
    DS[1]->SetLegend( prefix + "_H2" );
    DS[2]->SetLegend( prefix + "_H3" );
    SD[0]->SetLegend( "sd(" + prefix + "_H1)" );
    SD[1]->SetLegend( "sd(" + prefix + "_H2)" );
    SD[2]->SetLegend( "sd(" + prefix + "_H3)" );
    for (unsigned int nn = 0; nn != MAX_H_; nn++) {
      DS[nn]->SetDim(Dimension::X, Xaxis);
      SD[nn]->SetDim(Dimension::X, Xaxis);
    }

    // Loop over carbons in chain
    int pos = 0;
    for (ChainType::const_iterator it = Chains_[idx].begin(); it != Chains_[idx].end(); ++it)
    {
      if (debug_ > 0)
        mprintf("\t%s %s (%u)", resName, it->name(), it->Nvals());
      if (it->Nvals() > 0) {
        double avg, stdev;
        for (unsigned int i = 0; i != MAX_H_; i++) {
          //mprintf(" %15s", setup.Top().TruncResAtomName(site->Hidx(i)).c_str());
          if ( i < it->NumH() ) {
            if (report_p2_)
              avg = it->Avg(i, stdev);
            else
              avg = fabs(it->Avg(i, stdev));
          }
          if (debug_ > 0)
            mprintf("  %10.7f %10.7f  ", avg, stdev);
          DS[i]->Add(pos, &avg);
          SD[i]->Add(pos, &stdev);
        }
        if (debug_ > 0) mprintf("\n");
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
  double dnvals = (double)nvals_;
  double avg = sum_[i] / dnvals;
  stdev = sum2_[i] / dnvals;
  stdev -= (avg * avg);
  if (stdev > 0.0)
    stdev = sqrt( stdev );
  else
    stdev = 0.0;
  return avg;
}
