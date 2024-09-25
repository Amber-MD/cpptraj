#include "VDW_LongRange_Correction.h"
#include "../CpptrajStdio.h"
#include "../Topology.h"

using namespace Cpptraj::Energy;

VDW_LongRange_Correction::VDW_LongRange_Correction() :
  Vdw_Recip_term_(0),
  debug_(0)
{}

void VDW_LongRange_Correction::SetDebug(int debugIn) { debug_ = debugIn; }

/** Determine VDW long range correction prefactor. */
int VDW_LongRange_Correction::Setup_VDW_Correction(Topology const& topIn,
                                                    AtomMask const& maskIn)
{
  Vdw_Recip_term_ = 0.0;
  NonbondParmType const* NB_ = static_cast<NonbondParmType const*>( &(topIn.Nonbond()) );
  if (!NB_->HasNonbond()) {
    mprinterr("Error: '%s' has no nonbonded parameters. Cannot calculate VDW correction.\n",
            topIn.c_str());
    return 1;
  }
  // Count the number of each unique nonbonded type.
  N_vdw_type_.assign( NB_->Ntypes(), 0 );
  atype_vdw_recip_terms_.clear();
  atype_vdw_recip_terms_.reserve( N_vdw_type_.size() );
  vdw_type_.clear();
  vdw_type_.reserve( maskIn.Nselected() );
  for (AtomMask::const_iterator atm = maskIn.begin(); atm != maskIn.end(); ++atm)
  {
    N_vdw_type_[ topIn[*atm].TypeIndex() ]++;
    vdw_type_.push_back( topIn[*atm].TypeIndex() );
  }
  if (debug_ > 0) {
    mprintf("DEBUG: %zu VDW types.\n", N_vdw_type_.size());
    for (Iarray::const_iterator it = N_vdw_type_.begin(); it != N_vdw_type_.end(); ++it)
      mprintf("\tType %li = %i\n", it-N_vdw_type_.begin(), *it);
  }
  // Determine correction term from types and LJ B parameters
  for (unsigned int itype = 0; itype != N_vdw_type_.size(); itype++)
  {
    double atype_vdw_term = 0.0; // term for each nonbond atom type
    unsigned int offset = N_vdw_type_.size() * itype;
    for (unsigned int jtype = 0; jtype != N_vdw_type_.size(); jtype++)
    {
      unsigned int idx = offset + jtype;
      int nbidx = NB_->NBindex()[ idx ];
      if (nbidx > -1) {
        atype_vdw_term += N_vdw_type_[itype] * N_vdw_type_[jtype] * NB_->NBarray()[ nbidx ].B();

        Vdw_Recip_term_ += N_vdw_type_[itype] * N_vdw_type_[jtype] * NB_->NBarray()[ nbidx ].B();
      }
    }
    atype_vdw_recip_terms_.push_back(atype_vdw_term);
  }
  return 0;
}

/** Decomposed long range VDW correction */
double VDW_LongRange_Correction::Vdw_Decomp_Correction(std::vector<double>& atom_ecorr,
                                                       double cutoff_, double volume)
const
{
  atom_ecorr.resize( vdw_type_.size() );
  double prefac = Constants::TWOPI / (3.0*volume*cutoff_*cutoff_*cutoff_);
  double e_vdwr = -prefac * Vdw_Recip_term_;

  for (unsigned int i = 0; i != vdw_type_.size(); i++)
  {
    int v_type = vdw_type_[i];

    double term = atype_vdw_recip_terms_[v_type] / N_vdw_type_[v_type];

    atom_ecorr[i] = -prefac * term;
  }
  return e_vdwr;
}
