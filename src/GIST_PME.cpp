#ifdef LIBPME
#include <algorithm>
#include "GIST_PME.h"
#include "Frame.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include "ParameterTypes.h"
#include "Topology.h"

/** CONSTRUCTOR */
GIST_PME::GIST_PME() :
  NeighborCut2_(0)
{}

/** Strings corresponding to interaction type. */
const char* GIST_PME::InteractionTypeStr_[] = {
  "Other",
  "Solute0_Ongrid1",
  "Solvent0_Ongrid1",
  "Solute1_Ongrid0",
  "Solvent1_Ongrid0",
  "Both_Ongrid"
};

/** Setup up GIST PME calculation. Currently must be run on all atoms. */
int GIST_PME::Setup_PME_GIST(Topology const& topIn, unsigned int nthreads, double NeighborCut2_in)
{
  NeighborCut2_ = NeighborCut2_in;
  // Select everything
  allAtoms_ = AtomMask(0, topIn.Natom());
  // Set up PME
  if (Setup( topIn, allAtoms_ )) {
    mprinterr("Error: GIST PME setup failed.\n");
    return 1;
  }

  // Allocate voxel arrays
  int natoms = topIn.Natom(); // TODO unsigned int?
  E_vdw_direct_.resize( nthreads );
  E_elec_direct_.resize( nthreads );
  for (unsigned int t = 0; t != nthreads; t++)
  {
    E_vdw_direct_[t].assign(natoms, 0);
    E_elec_direct_[t].assign(natoms, 0);
  }
  // NOTE: Always allocate these arrays even if not using. This allows
  //       E_of_atom() function to work in all cases.
  //if (lw_coeff_ > 0) {
    E_vdw_self_.assign(natoms, 0);
    E_vdw_recip_.assign(natoms, 0);
  //} else
    E_vdw_lr_cor_.assign(natoms, 0);
  E_elec_self_.assign(natoms, 0);
  E_elec_recip_.assign(natoms, 0);

  e_potentialD_ = MatType(natoms, 4);
  e_potentialD_.setConstant(0.0);
  return 0;
}

/// \return Sum of elements in given array
static inline double SumDarray(std::vector<double> const& arr) {
  double sum = 0.0;
  for (std::vector<double>::const_iterator it = arr.begin(); it != arr.end(); ++it)
    sum += *it;
  return sum;
}


/** Calculate full nonbonded energy with PME Used for GIST, adding 6 arrays to store the
 * decomposed energy terms for every atom
 */
int GIST_PME::CalcNonbondEnergy_GIST(Frame const& frameIn,
                                     std::vector<int> const& atom_voxel,
                                     std::vector<bool> const& atomIsSolute,
                                     std::vector<bool> const& atomIsSolventO,
                                     std::vector<Darray>& E_UV_VDW_in,
                                     std::vector<Darray>& E_UV_Elec_in,
                                     std::vector<Darray>& E_VV_VDW_in,
                                     std::vector<Darray>& E_VV_Elec_in,
                                     std::vector<Farray>& Neighbor_in)
{
  t_total_.Start();
# ifdef DEBUG_GIST_PME
  mprintf("DEBUG: Beginning GIST PME nonbonded energy calculation for %zu atoms.\n", atom_voxel.size());
  for (unsigned int ii = 0; ii != atom_voxel.size(); ii++) {
    mprintf("\tAt %10u  voxel= %10i  isSolute= %1i\n", ii, atom_voxel[ii], (int)atomIsSolute[ii]);
  }
# endif
  // Elec. self
  double volume = frameIn.BoxCrd().CellVolume();
  std::fill(E_elec_self_.begin(), E_elec_self_.end(), 0);
  // decomposed E_elec_self for GIST
  double e_self = Self_GIST( volume, E_elec_self_, atom_voxel, atomIsSolute, E_VV_Elec_in[0] );

  // Create pairlist
  int retVal = pairList_.CreatePairList(frameIn, frameIn.BoxCrd().UnitCell(), frameIn.BoxCrd().FracCell(), allAtoms_);
  if (retVal != 0) {
    mprinterr("Error: GIST_PME: Grid setup failed.\n");
    return 1;
  }

  // TODO make more efficient
  coordsD_.clear();
  for (AtomMask::const_iterator atm = allAtoms_.begin(); atm != allAtoms_.end(); ++atm) {
    const double* XYZ = frameIn.XYZ( *atm );
    coordsD_.push_back( XYZ[0] );
    coordsD_.push_back( XYZ[1] );
    coordsD_.push_back( XYZ[2] );
  }

  unsigned int natoms = E_elec_self_.size();

  // Recip Potential for each atom
  e_potentialD_.setConstant(0.0);

  double e_recip = Recip_ParticleMesh_GIST( frameIn.BoxCrd(), atom_voxel, atomIsSolute,
                                            E_UV_Elec_in[0], E_VV_Elec_in[0] );

  // vdw potential for each atom 

  // TODO branch
  double e_vdw_lr_correction;
  double e_vdw6self, e_vdw6recip;
  if (lw_coeff_ > 0.0) {
    std::fill(E_vdw_self_.begin(), E_vdw_self_.end(), 0);
    e_vdw6self = Self6_GIST(E_vdw_self_);

    MatType vdw_potentialD(natoms, 4);
    vdw_potentialD.setConstant(0.0);
    e_vdw6recip = LJ_Recip_ParticleMesh_GIST( frameIn.BoxCrd(), vdw_potentialD );

    for(unsigned int j=0; j < natoms; j++){

      E_vdw_recip_[j] = 0.5 * (Cparam_[j] * vdw_potentialD(j,0)); // Split the energy by half

      // not sure vdw_recip for each atom can be calculated like this?
 
    } // by default, this block of code will not be excuted since the default lw_coeff_=-1


    if (debug_ > 0) {
      mprintf("DEBUG: gistpme e_vdw6self = %16.8f\n", e_vdw6self);
      mprintf("DEBUG: gistpme Evdwrecip = %16.8f\n", e_vdw6recip);
    }
    e_vdw_lr_correction = 0.0;
  } else {
    e_vdw6self = 0.0;
    e_vdw6recip = 0.0;
    std::fill(E_vdw_lr_cor_.begin(), E_vdw_lr_cor_.end(), 0);
    e_vdw_lr_correction = Vdw_Correction_GIST( volume, atom_voxel, atomIsSolute, E_UV_VDW_in[0], E_VV_VDW_in[0] );
    //mprintf("e_vdw_lr_correction: %f \n", e_vdw_lr_correction);
  }

  double e_vdw = 0.0;
  double e_direct = Direct_GIST( pairList_, e_vdw, atom_voxel, atomIsSolute, atomIsSolventO,
                                 E_UV_VDW_in, E_UV_Elec_in, E_VV_VDW_in, E_VV_Elec_in,
                                 Neighbor_in);

  //mprintf("e_elec_self: %f , e_elec_direct: %f, e_vdw6direct: %f \n", e_self, e_direct, e_vdw);

  if (debug_ > 0)
    mprintf("DEBUG: gistpme Eself= %20.10f   Erecip= %20.10f   Edirect= %20.10f  Evdw= %20.10f\n",
            e_self, e_recip, e_direct, e_vdw);

  //e_vdw += (e_vdw_lr_correction + e_vdw6self + e_vdw6recip);
  //e_elec = e_self + e_recip + e_direct;

# ifdef DEBUG_GIST_PME
  // DEBUG
  mprintf("DEBUG: gistpme voxel energies:\n");
  mprintf("\t%20s %20s %20s %20s\n", "E_UV_VDW", "E_UV_Elec", "E_VV_VDW", "E_VV_Elec");
  double sum_uv_vdw = 0;
  double sum_uv_elec = 0;
  double sum_vv_vdw = 0;
  double sum_vv_elec = 0;
  for (unsigned int vidx = 0; vidx != E_UV_VDW_in[0].size(); vidx++)
  {
    mprintf("\t%20.10g %20.10g %20.10g %20.10g\n",
            E_UV_VDW_in[0][vidx], E_UV_Elec_in[0][vidx], E_VV_VDW_in[0][vidx], E_VV_Elec_in[0][vidx]);
    sum_uv_vdw += E_UV_VDW_in[0][vidx];
    sum_uv_elec += E_UV_Elec_in[0][vidx];
    sum_vv_vdw += E_VV_VDW_in[0][vidx];
    sum_vv_elec += E_VV_Elec_in[0][vidx];
  }
  sum_uv_vdw /= 2.0;
  sum_uv_elec /= 2.0;
  sum_vv_vdw /= 2.0;
  sum_vv_elec /= 2.0;
  mprintf("DEBUG: Sums: UV_V= %20.10g  UV_E= %20.10g  VV_V= %20.10g  VV_E= %20.10g  Total= %20.10g\n",
           sum_uv_vdw, sum_uv_elec, sum_vv_vdw, sum_vv_elec,
           sum_uv_vdw + sum_uv_elec + sum_vv_vdw + sum_vv_elec);
# endif

  // DEBUG
  if (debug_ > 0) {
    // Calculate the sum of each terms
    double E_elec_direct_sum = SumDarray( E_elec_direct_[0] );
    double E_vdw_direct_sum = SumDarray( E_vdw_direct_[0] );

    double E_elec_self_sum   = SumDarray( E_elec_self_ );
    double E_elec_recip_sum  = SumDarray( E_elec_recip_ );
    mprintf("DEBUG: E_elec sums: self= %g  direct= %g  recip= %g\n", E_elec_self_sum, E_elec_direct_sum, E_elec_recip_sum);

    double E_vdw_self_sum   = SumDarray( E_vdw_self_ );
    double E_vdw_recip_sum  = SumDarray( E_vdw_recip_ );
    double E_vdw_lr_cor_sum = SumDarray( E_vdw_lr_cor_ );
    mprintf("DEBUG: E_vdw sums: self= %g  direct= %g  recip= %g  LR= %g\n", E_vdw_self_sum, E_vdw_direct_sum, E_vdw_recip_sum, E_vdw_lr_cor_sum);
  }

  t_total_.Stop();

  return 0;
}


/** Electrostatic self energy. This is the cancelling Gaussian plus the "neutralizing plasma". */
double GIST_PME::Self_GIST(double volume, Darray& atom_self,
                           std::vector<int> const& atom_voxel,
                           std::vector<bool> const& atomIsSolute,
                           Darray& e_vv_elec)
{
  t_self_.Start();
  double d0 = -ew_coeff_ * INVSQRTPI_;
  double ene = SumQ2() * d0;
//  mprintf("DEBUG: d0= %20.10f   ene= %20.10f\n", d0, ene);
  double factor = Constants::PI / (ew_coeff_ * ew_coeff_ * volume);
  double ee_plasma = -0.5 * factor * SumQ() * SumQ();

  for( unsigned int i=0; i< Charge_.size();i++)
  {
    // distribute the "neutrilizing plasma" to atoms equally
    atom_self[i]=Charge_[i]*Charge_[i]*d0 + ee_plasma/Charge_.size();
    if (atom_voxel[i] > -1 && !atomIsSolute[i]) {
      // Only do the self-terms for on-grid water; no exclusions for solute to on-grid solvent
      // NOTE: Multiply by 2 since other terms (direct/recip) are the 'full'
      //       terms (i.e. not divided between atoms) and the Action_GIST::CalcAvgVoxelEnergy()
      //       function will divide the voxels by a factor of 2.
      e_vv_elec[atom_voxel[i]] += (atom_self[i] * 2);
#     ifdef DEBUG_GIST_PME
      mprintf("DEBUG: gistpme vv self at %u eself= %20.10g\n", i,atom_self[i] );
#     endif
    }
  }

  ene += ee_plasma;
  t_self_.Stop();
  return ene;
}

/** Lennard-Jones self energy. for GIST */
double GIST_PME::Self6_GIST(Darray& atom_vdw_self) {
  t_self_.Start(); // TODO precalc
  double ew2 = lw_coeff_ * lw_coeff_;
  double ew6 = ew2 * ew2 * ew2;
  double c6sum = 0.0;
  for (Darray::const_iterator it = Cparam_.begin(); it != Cparam_.end(); ++it)
  {
    c6sum += ew6 * (*it * *it);

    atom_vdw_self.push_back(ew6 * (*it * *it)/12.0);

  }
  t_self_.Stop();
  return c6sum / 12.0;
}

/** PME recip calc for GIST to store decomposed recipical energy for every atom. */
double GIST_PME::Recip_ParticleMesh_GIST(Box const& boxIn,
                                         std::vector<int> const& atom_voxel,
                                         std::vector<bool> const& atomIsSolute,
                                         Darray& e_uv_elec, Darray& e_vv_elec)
{
  t_recip_.Start();
  // This essentially makes coordsD and chargesD point to arrays.
  MatType coordsD(&coordsD_[0], Charge_.size(), 3);
  MatType chargesD(&Charge_[0], Charge_.size(), 1);
  int nfft1 = Nfft(0);
  int nfft2 = Nfft(1);
  int nfft3 = Nfft(2);
  if ( DetermineNfft(nfft1, nfft2, nfft3, boxIn) ) {
    mprinterr("Error: Could not determine grid spacing.\n");
    return 0.0;
  }
  // Instantiate double precision PME object
  // Args: 1 = Exponent of the distance kernel: 1 for Coulomb
  //       2 = Kappa
  //       3 = Spline order
  //       4 = nfft1
  //       5 = nfft2
  //       6 = nfft3
  //       7 = scale factor to be applied to all computed energies and derivatives thereof
  //       8 = max # threads to use for each MPI instance; 0 = all available threads used.
  // NOTE: Scale factor for Charmm is 332.0716
  // NOTE: The electrostatic constant has been baked into the Charge_ array already.
  //auto pme_object = std::unique_ptr<PMEInstanceD>(new PMEInstanceD());
  pme_object_.setup(1, ew_coeff_, Order(), nfft1, nfft2, nfft3, 1.0, 0);
  // Sets the unit cell lattice vectors, with units consistent with those used to specify coordinates.
  // Args: 1 = the A lattice parameter in units consistent with the coordinates.
  //       2 = the B lattice parameter in units consistent with the coordinates.
  //       3 = the C lattice parameter in units consistent with the coordinates.
  //       4 = the alpha lattice parameter in degrees.
  //       5 = the beta lattice parameter in degrees.
  //       6 = the gamma lattice parameter in degrees.
  //       7 = lattice type
  pme_object_.setLatticeVectors(boxIn.Param(Box::X), boxIn.Param(Box::Y), boxIn.Param(Box::Z),
                                boxIn.Param(Box::ALPHA), boxIn.Param(Box::BETA), boxIn.Param(Box::GAMMA),
                                PMEInstanceD::LatticeType::XAligned);
  //double erecip = pme_object_.computeERec(0, chargesD, coordsD);
  double erecip = 0;
  pme_object_.computePRec(0,chargesD,coordsD,coordsD,1,e_potentialD_);
  for(unsigned int i =0; i < Charge_.size(); i++)
  {
    E_elec_recip_[i]=0.5 * Charge_[i] * e_potentialD_(i,0);
  }

  // For UV interaction, we need solute charges + positions and on-grid solvent positions.
  // For VV interaction, we need all solvent charges + positions and on-grid solvent positions.
  Darray Ucoords, Ucharges;
  Darray onGridCoords;
  std::vector<int> onGridAt;
  Darray Vcoords, Vcharges;
  unsigned int xidx = 0; // Index into coordsD_
  unsigned int n_on_grid = 0; // Number of solvent atoms on grid
  for (unsigned int atidx = 0; atidx != Charge_.size(); atidx++, xidx += 3)
  {
    if (atomIsSolute[atidx]) {
      Ucoords.push_back( coordsD_[xidx  ] );
      Ucoords.push_back( coordsD_[xidx+1] );
      Ucoords.push_back( coordsD_[xidx+2] );
      Ucharges.push_back( Charge_[atidx] );
    } else {
      Vcoords.push_back( coordsD_[xidx  ] );
      Vcoords.push_back( coordsD_[xidx+1] );
      Vcoords.push_back( coordsD_[xidx+2] );
      Vcharges.push_back( Charge_[atidx] );
      if (atom_voxel[atidx] > -1) {
        onGridCoords.push_back( coordsD_[xidx  ] );
        onGridCoords.push_back( coordsD_[xidx+1] );
        onGridCoords.push_back( coordsD_[xidx+2] );
        n_on_grid++;
        onGridAt.push_back( atidx );
      }
    }
  }
  MatType m_ux( &Ucoords[0],      Ucharges.size(),     3 );
  MatType m_uq( &Ucharges[0],     Ucharges.size(),     1 );
  MatType m_ox( &onGridCoords[0], n_on_grid,           3 );
  MatType m_vx( &Vcoords[0],      Vcharges.size(),     3 );
  MatType m_vq( &Vcharges[0],     Vcharges.size(),     1 );
  MatType ongrid_potentialD(n_on_grid, 4);
  // UV recip
  ongrid_potentialD.setZero();
  pme_object_.computePRec(0, m_uq, m_ux, m_ox, 1, ongrid_potentialD);
  for (unsigned int i = 0; i != onGridAt.size(); i++) {
#   ifdef DEBUG_GIST_PME
    mprintf("DEBUG: gistpme uv recip at %i erecip= %20.10g\n", onGridAt[i], Charge_[onGridAt[i]] * ongrid_potentialD(i, 0));
#   endif
    e_uv_elec[atom_voxel[onGridAt[i]]] += Charge_[onGridAt[i]] * ongrid_potentialD(i, 0);
  }
  // VV recip.
  // NOTE: If we do not zero the potential matrix, it will be summed into
  ongrid_potentialD.setZero();
  pme_object_.computePRec(0, m_vq, m_vx, m_ox, 1, ongrid_potentialD);
  for (unsigned int i = 0; i != onGridAt.size(); i++) {
#   ifdef DEBUG_GIST_PME
    mprintf("DEBUG: gistpme vv recip at %i erecip= %20.10g\n", onGridAt[i], Charge_[onGridAt[i]] * ongrid_potentialD(i, 0));
#   endif
    e_vv_elec[atom_voxel[onGridAt[i]]] += Charge_[onGridAt[i]] * ongrid_potentialD(i, 0);
  }

  t_recip_.Stop();
  return erecip;
}

/** The LJ PME reciprocal term for GIST*/ 
double GIST_PME::LJ_Recip_ParticleMesh_GIST(Box const& boxIn, MatType& potential)
{
  t_recip_.Start();
  int nfft1 = Nfft(0);
  int nfft2 = Nfft(1);
  int nfft3 = Nfft(2);
  if ( DetermineNfft(nfft1, nfft2, nfft3, boxIn) ) {
    mprinterr("Error: Could not determine grid spacing.\n");
    return 0.0;
  }

  MatType coordsD(&coordsD_[0], Charge_.size(), 3);
  MatType cparamD(&Cparam_[0], Cparam_.size(), 1);


  //auto pme_vdw = std::unique_ptr<PMEInstanceD>(new PMEInstanceD());
  pme_vdw_.setup(6, lw_coeff_, Order(), nfft1, nfft2, nfft3, -1.0, 0);
  pme_vdw_.setLatticeVectors(boxIn.Param(Box::X), boxIn.Param(Box::Y), boxIn.Param(Box::Z),
                             boxIn.Param(Box::ALPHA), boxIn.Param(Box::BETA), boxIn.Param(Box::GAMMA),
                             PMEInstanceD::LatticeType::XAligned);
  double evdwrecip = pme_vdw_.computeERec(0, cparamD, coordsD);
  pme_vdw_.computePRec(0,cparamD,coordsD,coordsD,1,potential);
  t_recip_.Stop();
  return evdwrecip;
}

/** Calculate full VDW long range correction from volume. */
double GIST_PME::Vdw_Correction_GIST(double volume,
                                     std::vector<int> const& atom_voxel,
                                     std::vector<bool> const& atomIsSolute,
                                     Darray& e_uv_vdw, Darray& e_vv_vdw)
{
  double prefac = Constants::TWOPI / (3.0*volume*cutoff_*cutoff_*cutoff_);
  //mprintf("VDW correction prefac: %.15f \n", prefac);
  double e_vdwr = -prefac * Vdw_Recip_Term();

  //mprintf("Cparam size: %i \n",Cparam_.size());

  //mprintf("volume of the unit cell: %f", volume);


  for ( unsigned int i = 0; i != vdw_type_.size(); i++)
  {
    double term(0);

    int v_type=vdw_type_[i];

    //v_type is the vdw_type of atom i, each atom has a atom type

    // atype_vdw_recip_terms_[vdw_type[i]] is the total vdw_recep_term for this v_type

    //N_vdw_type_[v_type] is the total atom number belongs to this v_type



    term = atype_vdw_recip_terms_[v_type] / N_vdw_type_[v_type];

    //mprintf("for i = %i,vdw_type = %i, Number of atoms in this vdw_type_= %i, Total vdw_recip_terms_ for this type= %f,  vdw_recip_term for atom i= %f \n",i,vdw_type_[i],N_vdw_type_[vdw_type_[i]],atype_vdw_recip_terms_[vdw_type_[i]],term);

    E_vdw_lr_cor_[i]= -prefac * term ;
    //mprintf("atom e_vdw_lr_cor: %f \n", -prefac* atom_vdw_recip_terms_[i]);
    int at_voxel = atom_voxel[i];
    if (at_voxel > -1) {
      // NOTE: Multiply by 2 since other terms (direct/recip) are the 'full'
      //       terms (i.e. not divided between atoms) and the Action_GIST::CalcAvgVoxelEnergy()
      //       function will divide the voxels by a factor of 2.
      if (atomIsSolute[i]) {
        e_uv_vdw[at_voxel] += (E_vdw_lr_cor_[i]*2);
#       ifdef DEBUG_GIST_PME
        mprintf("DEBUG: gistpme uv vdwLR at %u elr= %20.10g\n", i, E_vdw_lr_cor_[i]);
#       endif
      } else {
        e_vv_vdw[at_voxel] += (E_vdw_lr_cor_[i]*2);
#       ifdef DEBUG_GIST_PME
        mprintf("DEBUG: gistpme vv vdwLR at %u elr= %20.10g\n", i, E_vdw_lr_cor_[i]);
#       endif
      }
    }
  }

  if (debug_ > 0) mprintf("DEBUG: gistpme Vdw correction %20.10f\n", e_vdwr);
  return e_vdwr;
}

/** Direct space routine for GIST. */
double GIST_PME::Direct_GIST(PairList const& PL, double& evdw_out,
                             std::vector<int> const& atom_voxel,
                             std::vector<bool> const& atomIsSolute,
                             std::vector<bool> const& atomIsSolventO,
                             std::vector<Darray>& E_UV_VDW_in,
                             std::vector<Darray>& E_UV_Elec_in,
                             std::vector<Darray>& E_VV_VDW_in,
                             std::vector<Darray>& E_VV_Elec_in,
                             std::vector<Farray>& Neighbor_in)
            
{
  if (lw_coeff_ > 0.0)
    return Direct_VDW_LJPME_GIST(PL, evdw_out);
  else
    return Direct_VDW_LongRangeCorrection_GIST(PL, evdw_out, atom_voxel, atomIsSolute, atomIsSolventO,
                                               E_UV_VDW_in, E_UV_Elec_in,
                                               E_VV_VDW_in, E_VV_Elec_in,
                                               Neighbor_in);
}

/** Set interaction type and optionally on-grid solvent voxel index.
  * Interaction type: 0=other, 1=solute to on-grid solvent, 2=any solvent to on-grid solvent.
  */
static inline void determineInteractionType(int& interactionType, int& onGridVoxelIdx,
                                            int voxel0, bool isSolute0,
                                            int voxel1, bool isSolute1)
{
  interactionType = 0;
  onGridVoxelIdx = -1;
  bool onGridSolvent0 = (voxel0 > -1 && !isSolute0);
  bool onGridSolvent1 = (voxel1 > -1 && !isSolute1);
  // We only care if at least 1 atom corresponds to an on-grid solvent molecule
  if (onGridSolvent1) {
    onGridVoxelIdx = voxel1;
    if (isSolute0)
      // Solute0 to on-grid solvent 1
      interactionType = 1;
    else
      // Solvent0 to on-grid solvent 1
      interactionType = 2;
  } else if (onGridSolvent0) {
    onGridVoxelIdx = voxel0;
    if (isSolute1)
      // Solute1 to on-grid solvent 0
      interactionType = 1;
    else
      // Solvent1 to on-grid solvent 0
      interactionType = 2;
  }
}

/** \return Interaction type between atoms given grid voxels and solute status. */
GIST_PME::InteractionType GIST_PME::determineInteractionType(int voxel0, bool isSolute0,
                                                             int voxel1, bool isSolute1)
{
  InteractionType interactionType = OTHER;
  bool onGridSolvent0 = (voxel0 > -1 && !isSolute0);
  bool onGridSolvent1 = (voxel1 > -1 && !isSolute1);
  
  if (onGridSolvent1) {
    if (onGridSolvent0)
      // Both are on the grid
      interactionType = BOTH_ONGRID;
    else if (isSolute0)
      // Solute0 to on-grid solvent 1
      interactionType = SOLUTE0_ONGRID1;
    else
      // Solvent0 to on-grid solvent 1
      interactionType = SOLVENT0_ONGRID1;
  } else if (onGridSolvent0) {
    if (onGridSolvent1)
      // Both are on the grid
      interactionType = BOTH_ONGRID;
    else if (isSolute1)
      // Solute1 to on-grid solvent 0
      interactionType = SOLUTE1_ONGRID0;
    else
      // Solvent1 to on-grid solvent 0
      interactionType = SOLVENT1_ONGRID0;
  }
  return interactionType;
}

/** Nonbond energy kernel.
  * \param Eelec Total electrostatic energy; will be incremented.
  * \param Evdw Total van der Waals energy; will be incremented.
  * \param rij2 Distance between the atoms squared.
  * \param idx0 Index of the first atom.
  * \param idx1 Index of the second atom.
  * \param e_elec_direct Direct space electrostatic energy array (natoms).
  * \param e_vdw_direct Direct space van der Waals energy array (natoms).
  * \param interactionType 0=other, 1=solute to on-grid solvent, 2=any solvent to on-grid solvent.
  * \param onGridVoxelIdx Index of on-grid solvent for e_??_X arrays for interactionType > 0.
  * \param e_uv_vdw Solute to on-grid solvent van der Waals array (nvoxels).
  * \param e_uv_elec Solute to on-grid solvent electrostatics array (nvoxels).
  * \param e_vv_vdw Any solvent to on-grid solvent var der Waals array (nvoxels).
  * \param e_vv_elec Any solvent to on-grid solvent electrostatics array (nvoxels).
  */
void GIST_PME::Ekernel_NB(double& Eelec, double& Evdw,
                          double rij2,
                          double q0, double q1,
                          int idx0, int idx1,
                          double* e_elec_direct, double* e_vdw_direct,
                          InteractionType interactionType, int voxel0, int voxel1,
                          double* e_uv_vdw, double* e_uv_elec,
                          double* e_vv_vdw, double* e_vv_elec)
//const Cannot be const because of the timer
{
                double rij = sqrt( rij2 );
                double qiqj = q0 * q1;
#               ifndef _OPENMP
                t_erfc_.Start();
#               endif
                //double erfc = erfc_func(ew_coeff_ * rij);
                double erfc = ErfcFxn(ew_coeff_ * rij);
#               ifndef _OPENMP
                t_erfc_.Stop();
#               endif
                double e_elec = qiqj * erfc / rij;
                Eelec += e_elec;

                //add the e_elec to E_elec_direct[it0->Idx] and E_elec_direct[it1-> Idx], 0.5 to avoid double counting
                 
                e_elec_direct[idx0] += 0.5 * e_elec;
                e_elec_direct[idx1] += 0.5 * e_elec;

                // DEBUG
#               ifdef DEBUG_GIST_PME
                mprintf("DEBUG: gistpme direct itype='%s' e_elec= %20.10f\n", InteractionTypeStr_[interactionType], e_elec);
#               endif
                if (interactionType == SOLUTE0_ONGRID1)
                  e_uv_elec[voxel1] += e_elec;
                else if (interactionType == SOLVENT0_ONGRID1)
                  e_vv_elec[voxel1] += e_elec;
                else if (interactionType == SOLUTE1_ONGRID0)
                  e_uv_elec[voxel0] += e_elec;
                else if (interactionType == SOLVENT1_ONGRID0)
                  e_vv_elec[voxel0] += e_elec;
                else if (interactionType == BOTH_ONGRID) {
                  e_vv_elec[voxel0] += e_elec;
                  e_vv_elec[voxel1] += e_elec;
                }

                int nbindex = NB().GetLJindex(TypeIdx(idx0),
                                              TypeIdx(idx1));
                if (nbindex > -1) {
                  double vswitch = SwitchFxn(rij2, cut2_0_, cut2_);
                  NonbondType const& LJ = NB().NBarray()[ nbindex ];
                  double r2    = 1.0 / rij2;
                  double r6    = r2 * r2 * r2;
                  double r12   = r6 * r6;
                  double f12   = LJ.A() * r12;  // A/r^12
                  double f6    = LJ.B() * r6;   // B/r^6
                  double e_vdw = f12 - f6;      // (A/r^12)-(B/r^6)
                  double e_vdw_term = e_vdw * vswitch;
                  Evdw += e_vdw_term;

                  // add the vdw energy to e_vdw_direct, divide the energy evenly to each atom

                  e_vdw_direct[idx0] += 0.5 * e_vdw_term;
                  e_vdw_direct[idx1] += 0.5 * e_vdw_term;
#                 ifdef DEBUG_GIST_PME
                  mprintf("DEBUG: gistpme direct itype='%s' e_vdw_term= %20.10f\n", InteractionTypeStr_[interactionType], e_vdw_term);
#                 endif
                  if (interactionType == SOLUTE0_ONGRID1)
                    e_uv_vdw[voxel1] += e_vdw_term;
                  else if (interactionType == SOLVENT0_ONGRID1)
                    e_vv_vdw[voxel1] += e_vdw_term;
                  else if (interactionType == SOLUTE1_ONGRID0)
                    e_uv_vdw[voxel0] += e_vdw_term;
                  else if (interactionType == SOLVENT1_ONGRID0)
                    e_vv_vdw[voxel0] += e_vdw_term;
                  else if (interactionType == BOTH_ONGRID) {
                    e_vv_vdw[voxel0] += e_vdw_term;
                    e_vv_vdw[voxel1] += e_vdw_term;
                  }

                  //mprintf("PVDW %8i%8i%20.6f%20.6f\n", ta0+1, ta1+1, e_vdw, r2);
#                 ifdef CPPTRAJ_EKERNEL_LJPME
                  // LJ PME direct space correction
                  double kr2 = lw_coeff_ * lw_coeff_ * rij2;
                  double kr4 = kr2 * kr2;
                  //double kr6 = kr2 * kr4;
                  double expterm = exp(-kr2);
                  double Cij = Cparam_[idx0] * Cparam_[idx1];
                  Eljpme_correction += (1.0 - (1.0 +  kr2 + kr4/2.0)*expterm) * r6 * vswitch * Cij;
#                 endif
                }
}

/** Adjust energy kernel. */
void GIST_PME::Ekernel_Adjust(double& e_adjust,
                              double rij2,
                              double q0, double q1,
                              int idx0, int idx1,
                              double* e_elec_direct,
                              InteractionType interactionType, int voxel0, int voxel1,
                              double* e_uv_elec, double* e_vv_elec)
{

              double adjust = AdjustFxn(q0,q1,sqrt(rij2));

              e_adjust += adjust;

              //add the e_elec to E_elec_direct[it0->Idx] and E_elec_direct[it1-> Idx], 0.5 to avoid double counting

              e_elec_direct[idx0] += 0.5 * adjust;
              e_elec_direct[idx1] += 0.5 * adjust;
#             ifdef DEBUG_GIST_PME
              mprintf("DEBUG: gistpme adjust itype='%s' adjust= %20.10f\n", InteractionTypeStr_[interactionType], adjust);
#             endif
              if (interactionType == SOLUTE0_ONGRID1)
                e_uv_elec[voxel1] += adjust;
              else if (interactionType == SOLVENT0_ONGRID1)
                e_vv_elec[voxel1] += adjust;
              else if (interactionType == SOLUTE1_ONGRID0)
                e_uv_elec[voxel0] += adjust;
              else if (interactionType == SOLVENT1_ONGRID0)
                e_vv_elec[voxel0] += adjust;
              else if (interactionType == BOTH_ONGRID) {
                e_vv_elec[voxel0] += adjust;
                e_vv_elec[voxel1] += adjust;
              }

#             ifdef CPPTRAJ_EKERNEL_LJPME
              // LJ PME direct space correction
              // NOTE: Assuming excluded pair is within cutoff
              double kr2 = lw_coeff_ * lw_coeff_ * rij2;
              double kr4 = kr2 * kr2;
              //double kr6 = kr2 * kr4;
              double expterm = exp(-kr2);
              double r4 = rij2 * rij2;
              double r6 = rij2 * r4;
              double Cij = Cparam_[it0->Idx()] * Cparam_[it1->Idx()];
              Eljpme_correction_excl += (1.0 - (1.0 +  kr2 + kr4/2.0)*expterm) / r6 * Cij;
#             endif
}

static inline void neighborCalc(float* neighbor, double rij2, double NeighborCut2, int a1_voxel, int a2_voxel)
{
  if (rij2 < NeighborCut2) {
    if (a1_voxel > -1) neighbor[a1_voxel] += 1.0;
    if (a2_voxel > -1) neighbor[a2_voxel] += 1.0;
  }
}

/** Direct space calculation with long range VDW correction for GIST. */
double GIST_PME::Direct_VDW_LongRangeCorrection_GIST(PairList const& PL, double& evdw_out,
                                                     std::vector<int> const& atom_voxel,
                                                     std::vector<bool> const& atomIsSolute,
                                                     std::vector<bool> const& atomIsSolventO,
                                                     std::vector<Darray>& E_UV_VDW_in,
                                                     std::vector<Darray>& E_UV_Elec_in,
                                                     std::vector<Darray>& E_VV_VDW_in,
                                                     std::vector<Darray>& E_VV_Elec_in,
                                                     std::vector<Farray>& Neighbor_in)
{
  // e_vdw_direct only count the interaction water molecule involved( sw, ww)
  // e_vdw_all_direct, count all interactions( sw, ww, ss)

  t_direct_.Start();
  double Eelec = 0.0;
  double e_adjust = 0.0;
  double Evdw = 0.0;
  int cidx;

  //mprintf("Entering Direct_VDW_LongrangeCorrection_GIST function \n");
# ifndef _OPENMP
  // Zero when not openmp
  std::fill(E_vdw_direct_[0].begin(), E_vdw_direct_[0].end(), 0);
  std::fill(E_elec_direct_[0].begin(), E_elec_direct_[0].end(), 0);
# endif
  double* e_vdw_direct  = &(E_vdw_direct_[0][0]);
  double* e_elec_direct = &(E_elec_direct_[0][0]);
  double* e_uv_vdw  = &(E_UV_VDW_in[0][0]);
  double* e_uv_elec = &(E_UV_Elec_in[0][0]);
  double* e_vv_vdw  = &(E_VV_VDW_in[0][0]);
  double* e_vv_elec = &(E_VV_Elec_in[0][0]);
  float*  neighbor  = &(Neighbor_in[0][0]);
  // Pair list loop
# ifdef _OPENMP
  int mythread;
# pragma omp parallel private(cidx, mythread, e_vdw_direct, e_elec_direct, e_uv_vdw, e_uv_elec, e_vv_vdw, e_vv_elec, neighbor) reduction(+: Eelec, Evdw, e_adjust)
  {
  mythread = omp_get_thread_num();
  std::fill(E_vdw_direct_[mythread].begin(), E_vdw_direct_[mythread].end(), 0);
  std::fill(E_elec_direct_[mythread].begin(), E_elec_direct_[mythread].end(), 0);
  e_vdw_direct  = &(E_vdw_direct_[mythread][0]);
  e_elec_direct = &(E_elec_direct_[mythread][0]);
  e_uv_vdw = &(E_UV_VDW_in[mythread][0]);
  e_uv_elec = &(E_UV_Elec_in[mythread][0]);
  e_vv_vdw = &(E_VV_VDW_in[mythread][0]);
  e_vv_elec = &(E_VV_Elec_in[mythread][0]);
  neighbor = &(Neighbor_in[mythread][0]);
# pragma omp for
# endif
//# include "PairListLoop.h"
  for (cidx = 0; cidx < PL.NGridMax(); cidx++)
  {
    PairList::CellType const& thisCell = PL.Cell( cidx );
    if (thisCell.NatomsInGrid() > 0)
    {
      // cellList contains this cell index and all neighbors.
      PairList::Iarray const& cellList = thisCell.CellList();
      // transList contains index to translation for the neighbor.
      PairList::Iarray const& transList = thisCell.TransList();
      // Loop over all atoms of thisCell.
      for (PairList::CellType::const_iterator it0 = thisCell.begin();
                                              it0 != thisCell.end(); ++it0)
      {
        Vec3 const& xyz0 = it0->ImageCoords();
        double q0 = Charge_[it0->Idx()];
        // The voxel # of it0
        int it0_voxel = atom_voxel[it0->Idx()];
        bool it0_solute = atomIsSolute[it0->Idx()];
        bool it0_solventO = atomIsSolventO[it0->Idx()];

#       ifdef DEBUG_PAIRLIST
        mprintf("DBG: Cell %6i (%6i atoms):\n", cidx+1, thisCell.NatomsInGrid());
#       endif
        // Exclusion list for this atom
        ExclusionArray::ExListType const& excluded = Excluded()[it0->Idx()];
        // Calc interaction of atom to all other atoms in thisCell.
        for (PairList::CellType::const_iterator it1 = it0 + 1;
                                                it1 != thisCell.end(); ++it1)
        {
          int it1_voxel = atom_voxel[it1->Idx()];
          bool it1_solute = atomIsSolute[it1->Idx()];
          bool it1_solventO = atomIsSolventO[it1->Idx()];
          if (it0_voxel > -1 || it1_voxel > -1) {
            Vec3 const& xyz1 = it1->ImageCoords();
            double q1 = Charge_[it1->Idx()];
            Vec3 dxyz = xyz1 - xyz0;
            double rij2 = dxyz.Magnitude2();
#           ifdef DEBUG_PAIRLIST
            mprintf("\tAtom %6i to atom %6i (%f)\n", it0->Idx()+1, it1->Idx()+1, sqrt(rij2));
#           endif
            InteractionType interactionType = determineInteractionType(it0_voxel, it0_solute, it1_voxel, it1_solute);
            // If atom excluded, calc adjustment, otherwise calc elec. energy.
            if (excluded.find( it1->Idx() ) == excluded.end())
            {
              if ( rij2 < cut2_ ) {
                Ekernel_NB(Eelec, Evdw, rij2, q0, q1, it0->Idx(), it1->Idx(), e_elec_direct, e_vdw_direct,
                             interactionType, it0_voxel, it1_voxel, 
                             e_uv_vdw, e_uv_elec, e_vv_vdw, e_vv_elec);
              }
            } else {
                Ekernel_Adjust(e_adjust, rij2, q0, q1, it0->Idx(), it1->Idx(), e_elec_direct,
                               interactionType, it0_voxel, it1_voxel, e_uv_elec, e_vv_elec);
            }
            if (it0_solventO && it1_solventO) {
              neighborCalc(neighbor, rij2, NeighborCut2_, it0_voxel, it1_voxel);
            }
          } // END at least 1 voxel on grid
        } // END loop over other atoms in thisCell
        // Loop over all neighbor cells
        for (unsigned int nidx = 1; nidx != cellList.size(); nidx++)
        {
          PairList::CellType const& nbrCell = PL.Cell( cellList[nidx] );
#         ifdef DEBUG_PAIRLIST
          if (nbrCell.NatomsInGrid()>0) mprintf("\tto neighbor cell %6i\n", cellList[nidx]+1);
#         endif
          // Translate vector for neighbor cell
          Vec3 const& tVec = PL.TransVec( transList[nidx] );
          //mprintf("\tNEIGHBOR %i (idxs %i - %i)\n", nbrCell, beg1, end1);
          // Loop over every atom in nbrCell
          for (PairList::CellType::const_iterator it1 = nbrCell.begin();
                                                  it1 != nbrCell.end(); ++it1)
          {
            int it1_voxel = atom_voxel[it1->Idx()];
            bool it1_solute = atomIsSolute[it1->Idx()];
            bool it1_solventO = atomIsSolventO[it1->Idx()];
            if (it0_voxel > -1 || it1_voxel > -1) {
              Vec3 const& xyz1 = it1->ImageCoords();
              double q1 = Charge_[it1->Idx()];
              Vec3 dxyz = xyz1 + tVec - xyz0;
              double rij2 = dxyz.Magnitude2();
#             ifdef DEBUG_PAIRLIST
              mprintf("\t\tAtom %6i to atom %6i (%f)\n", it0->Idx()+1, it1->Idx()+1, sqrt(rij2));
#             endif
              //mprintf("\t\tNbrAtom %06i\n",atnum1);
              InteractionType interactionType = determineInteractionType(it0_voxel, it0_solute, it1_voxel, it1_solute);
              // If atom excluded, calc adjustment, otherwise calc elec. energy.
              // TODO Is there better way of checking this?
              if (excluded.find( it1->Idx() ) == excluded.end())
              {
                //mprintf("\t\t\tdist= %f\n", sqrt(rij2));
                if ( rij2 < cut2_ ) {
                  Ekernel_NB(Eelec, Evdw, rij2, q0, q1, it0->Idx(), it1->Idx(), e_elec_direct, e_vdw_direct,
                             interactionType, it0_voxel, it1_voxel, 
                             e_uv_vdw, e_uv_elec, e_vv_vdw, e_vv_elec);
                }
              } else {
                Ekernel_Adjust(e_adjust, rij2, q0, q1, it0->Idx(), it1->Idx(), e_elec_direct,
                               interactionType, it0_voxel, it1_voxel, e_uv_elec, e_vv_elec);
              }
              if (it0_solventO && it1_solventO) {
                neighborCalc(neighbor, rij2, NeighborCut2_, it0_voxel, it1_voxel);
              }
            } // END at least 1 voxel on the grid
          } // END loop over neighbor cell atoms
        
        } // END Loop over neighbor cells
        
      } // Loop over thisCell atoms

    } // END if thisCell is not empty
  } // Loop over cells

# ifdef _OPENMP
  } // END pragma omp parallel
  // Sum up contributions from individual threads into array 0
  for (unsigned int t = 1; t < E_vdw_direct_.size(); t++) {
    for (unsigned int i = 0; i != E_vdw_direct_[t].size(); i++) {
      E_vdw_direct_[0][i] += E_vdw_direct_[t][i];
      E_elec_direct_[0][i] += E_elec_direct_[t][i];
    }
  }
# endif
  t_direct_.Stop();
# ifdef DEBUG_PAIRLIST
  mprintf("DEBUG: gpme Elec                             = %16.8f\n", Eelec);
  mprintf("DEBUG: gpme Eadjust                          = %16.8f\n", e_adjust);
  mprintf("DEBUG: gpme LJ vdw                           = %16.8f\n", Evdw);
# endif
  evdw_out = Evdw;
  return Eelec + e_adjust;
}

/** Direct space calculation with LJ PME for GIST. */ // TODO enable
double GIST_PME::Direct_VDW_LJPME_GIST(PairList const& PL, double& evdw_out)
{
  mprinterr("Error: LJPME does not yet work with GIST.\n");
  return 0;
/*
  t_direct_.Start();
  double Eelec = 0.0;
  double e_adjust = 0.0;
  double Evdw = 0.0;
  double Eljpme_correction = 0.0;
  double Eljpme_correction_excl = 0.0;
  int cidx;
# define CPPTRAJ_EKERNEL_LJPME
# ifdef _OPENMP
# pragma omp parallel private(cidx) reduction(+: Eelec, Evdw, e_adjust, Eljpme_correction,Eljpme_correction_excl)
  {
# pragma omp for
# endif
//# include "PairListLoop.h"
  double Evdw_temp(0),Eljpme_correction_temp(0), Eljpme_correction_excl_temp(0);
  double Eelec_temp(0),E_adjust_temp(0);
  
  
  for (cidx = 0; cidx < PL.NGridMax(); cidx++)
  {
    PairList::CellType const& thisCell = PL.Cell( cidx );
    if (thisCell.NatomsInGrid() > 0)
    {
      // cellList contains this cell index and all neighbors.
      PairList::Iarray const& cellList = thisCell.CellList();
      // transList contains index to translation for the neighbor.
      PairList::Iarray const& transList = thisCell.TransList();
      // Loop over all atoms of thisCell.
      for (PairList::CellType::const_iterator it0 = thisCell.begin();
                                              it0 != thisCell.end(); ++it0)
      {
        Vec3 const& xyz0 = it0->ImageCoords();
        double q0 = Charge_[it0->Idx()];
#       ifdef DEBUG_PAIRLIST
        mprintf("DBG: Cell %6i (%6i atoms):\n", cidx+1, thisCell.NatomsInGrid());
#       endif
        // Exclusion list for this atom
        ExclusionArray::ExListType const& excluded = Excluded()[it0->Idx()];
        // Calc interaction of atom to all other atoms in thisCell.
        for (PairList::CellType::const_iterator it1 = it0 + 1;
                                                it1 != thisCell.end(); ++it1)
        {
          Vec3 const& xyz1 = it1->ImageCoords();
          double q1 = Charge_[it1->Idx()];
          Vec3 dxyz = xyz1 - xyz0;
          double rij2 = dxyz.Magnitude2();
#         ifdef DEBUG_PAIRLIST
          mprintf("\tAtom %6i to atom %6i (%f)\n", it0->Idx()+1, it1->Idx()+1, sqrt(rij2));
#         endif
          // If atom excluded, calc adjustment, otherwise calc elec. energy.
          if (excluded.find( it1->Idx() ) == excluded.end())
          {
            if ( rij2 < cut2_ ) {
#             include "EnergyKernel_Nonbond.h"
            }
          } else {
#           include "EnergyKernel_Adjust.h"
          }
        } // END loop over other atoms in thisCell
        // Loop over all neighbor cells
        for (unsigned int nidx = 1; nidx != cellList.size(); nidx++)
        {
          PairList::CellType const& nbrCell = PL.Cell( cellList[nidx] );
#         ifdef DEBUG_PAIRLIST
          if (nbrCell.NatomsInGrid()>0) mprintf("\tto neighbor cell %6i\n", cellList[nidx]+1);
#         endif
          // Translate vector for neighbor cell
          Vec3 const& tVec = PL.TransVec( transList[nidx] );
          //mprintf("\tNEIGHBOR %i (idxs %i - %i)\n", nbrCell, beg1, end1);
          // Loop over every atom in nbrCell
          for (PairList::CellType::const_iterator it1 = nbrCell.begin();
                                                  it1 != nbrCell.end(); ++it1)
          {
            Vec3 const& xyz1 = it1->ImageCoords();
            double q1 = Charge_[it1->Idx()];
            Vec3 dxyz = xyz1 + tVec - xyz0;
            double rij2 = dxyz.Magnitude2();
#           ifdef DEBUG_PAIRLIST
            mprintf("\t\tAtom %6i to atom %6i (%f)\n", it0->Idx()+1, it1->Idx()+1, sqrt(rij2));
#           endif
            //mprintf("\t\tNbrAtom %06i\n",atnum1);
            // If atom excluded, calc adjustment, otherwise calc elec. energy.
            // TODO Is there better way of checking this?
            if (excluded.find( it1->Idx() ) == excluded.end())
            {
              //mprintf("\t\t\tdist= %f\n", sqrt(rij2));
              if ( rij2 < cut2_ ) {
#               include "EnergyKernel_Nonbond.h"
              }
            } else {
#             include "EnergyKernel_Adjust.h"
            }
          } // END loop over neighbor cell atoms
        
        } // END Loop over neighbor cells
        

        e_vdw_direct[it0->Idx()]=(Evdw + Eljpme_correction + Eljpme_correction_excl)-(Evdw_temp + Eljpme_correction_temp + Eljpme_correction_excl_temp);
        
        e_elec_direct[it0->Idx()]=(Eelec + e_adjust) -(Eelec_temp + E_adjust_temp);

        Evdw_temp=Evdw;
        Eljpme_correction_temp=Eljpme_correction;
        Eljpme_correction_excl_temp=Eljpme_correction_excl;
        Eelec_temp=Eelec;
        E_adjust_temp=e_adjust;
      } // Loop over thisCell atoms
      


    } // END if thisCell is not empty
  } // Loop over cells

# ifdef _OPENMP
  } // END pragma omp parallel
# endif
# undef CPPTRAJ_EKERNEL_LJPME
  t_direct_.Stop();
# ifdef DEBUG_PAIRLIST
  mprintf("DEBUG: Elec                             = %16.8f\n", Eelec);
  mprintf("DEBUG: Eadjust                          = %16.8f\n", e_adjust);
  mprintf("DEBUG: LJ vdw                           = %16.8f\n", Evdw);
  mprintf("DEBUG: LJ vdw PME correction            = %16.8f\n", Eljpme_correction);
  mprintf("DEBUG: LJ vdw PME correction (excluded) = %16.8f\n", Eljpme_correction_excl);
# endif
  evdw_out = Evdw + Eljpme_correction + Eljpme_correction_excl;
  return Eelec + e_adjust;
*/
}
#endif /* LIBPME */
