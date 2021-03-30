#ifdef LIBPME
#include "GIST_PME.h"
#include "Frame.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include "ParameterTypes.h"

/** CONSTRUCTOR */
GIST_PME::GIST_PME() {}

/** Allocate internal arrays. */
int GIST_PME::AllocateArrays(unsigned int natoms, unsigned int nvoxels, unsigned int nthreads)
{
  // Allocate voxel arrays
  E_vdw_direct_.resize( nthreads );
  E_elec_direct_.resize( nthreads );
  for (unsigned int t = 0; t != nthreads; t++)
  {
    E_vdw_direct_[t].assign(nvoxels, 0);
    E_elec_direct_[t].assign(nvoxels, 0);
  }
  if (lw_coeff_ > 0) {
    E_vdw_self_.assign(nvoxels, 0);
    E_vdw_recip_.assign(nvoxels, 0);
  } else
    E_vdw_lr_cor_.assign(nvoxels, 0);
  E_elec_self_.assign(nvoxels, 0);
  E_elec_recip_.assign(nvoxels, 0);

  e_potentialD_ = MatType(natoms, 4);
  e_potentialD_.setConstant(0.0);
  return 0;
}

/** Calculate full nonbonded energy with PME Used for GIST, adding 6 arrays to store the
 * decomposed energy terms for every atom
 */
int GIST_PME::CalcNonbondEnergy_GIST(Frame const& frameIn, AtomMask const& maskIn,
                                      double& e_elec, double& e_vdw,
                                      Darray& e_vdw_direct,
                                      Darray& e_vdw_self,
                                      Darray& e_vdw_recip,
                                      Darray& e_vdw_lr_cor,
                                      Darray& e_elec_self,
                                      Darray& e_elec_direct,
                                      Darray& e_elec_recip,
                                      Iarray const& atom_voxel )
{
  t_total_.Start();
  double volume = frameIn.BoxCrd().CellVolume();

  //auto step0 = std::chrono::system_clock::now();


  double e_self = Self_GIST( volume, e_elec_self ); // decomposed E_elec_self for GIST

  //auto step1 = std::chrono::system_clock::now();

  //std::chrono::duration<double> d1 = step1 -step0;

  //mprintf("Eelec_self takes: %f seconds\n", d1.count());



  int retVal = pairList_.CreatePairList(frameIn, frameIn.BoxCrd().UnitCell(), frameIn.BoxCrd().FracCell(), maskIn);
  if (retVal != 0) {
    mprinterr("Error: Grid setup failed.\n");
    return 1;
  }

  // TODO make more efficient
  int idx = 0;
  coordsD_.clear();
  

  for (AtomMask::const_iterator atm = maskIn.begin(); atm != maskIn.end(); ++atm, ++idx) {
    const double* XYZ = frameIn.XYZ( *atm );
    coordsD_.push_back( XYZ[0] );
    coordsD_.push_back( XYZ[1] );
    coordsD_.push_back( XYZ[2] );
    
  }

  const int atom_num = coordsD_.size()/3;

  //mprintf("atom_num is %i \n", atom_num);

 //Potential for each atom
 //helpme::Matrix<double> e_potentialD[atom_num]={0};
  MatType e_potentialD(atom_num,4);

  e_potentialD.setConstant(0.0);
  //mprintf("The vaule at row 10,col 1: %f \n", e_potentialD(10,1));

  //Mat chargesD(&Charge_[0], Charge_.size(), 1);

  //Mat e_potentialD = chargesD.clone()
  

//  MapCoords(frameIn, ucell, recip, maskIn);
  double e_recip = Recip_ParticleMesh_GIST( frameIn.BoxCrd(), e_potentialD );

  //auto step2 = std::chrono::system_clock::now();


  //mprintf("e_recip finished, the e_recip value: %f \n", e_recip);

  //std::chrono::duration<double> d2 = step2 - step1;

  //mprintf("Eelec_recip takes: %f seconds \n", d2.count());


  //Darray* iterator it;
  //int atom_num =0;

  //idx=0;

  /**
  
  for ( Ewald::Darray *it=Charge_; ++it,++idx)
  {
    e_elec_recip.push_back(*it * e_potentialD(idx,1));
    //atom_num=atom_num+1;

  }
  **/

  for(int i =0; i < atom_num;i++)
  {
    e_elec_recip[i]=0.5 * Charge_[i] * e_potentialD(i,0);
  }

  //auto step3= std::chrono::system_clock::now();

  //std::chrono::duration<double> d3 = step3 -step2;
  //mprintf("decompose Eelec_recip takes: %f seconds \n", d3.count());


  

  // vdw potential for each atom 


  // TODO branch
  double e_vdw_lr_correction;
  double e_vdw6self, e_vdw6recip;
  if (lw_coeff_ > 0.0) {
    MatType vdw_potentialD(atom_num,4);

    e_vdw6self = Self6_GIST(e_vdw_self);

    e_vdw6recip = LJ_Recip_ParticleMesh_GIST( frameIn.BoxCrd(),vdw_potentialD );

    mprintf(" e_vdw6self: %f, e_vdw6recip: %f \n", e_vdw6self, e_vdw6recip);

    for(int j=0; j< atom_num;j++){

    e_vdw_recip[j]=0.5 * (Cparam_[j] * vdw_potentialD(j,0)); // Split the energy by half

     // not sure vdw_recip for each atom can be calculated like this?
 
    } // by default, this block of code will not be excuted since the default lw_coeff_=-1


    if (debug_ > 0) {
      mprintf("DEBUG: e_vdw6self = %16.8f\n", e_vdw6self);
      mprintf("DEBUG: Evdwrecip = %16.8f\n", e_vdw6recip);
    }
    e_vdw_lr_correction = 0.0;
  } else {
    e_vdw6self = 0.0;
    e_vdw6recip = 0.0;
    e_vdw_lr_correction = Vdw_Correction_GIST( volume, e_vdw_lr_cor );
    //mprintf("e_vdw_lr_correction: %f \n", e_vdw_lr_correction);
  }

  //auto step4= std::chrono::system_clock::now();

  //mprintf("vdw long range correction takes: %f seconds \n", (step4-step3).count());

  e_vdw = 0.0;
  double e_direct = Direct_GIST( pairList_, e_vdw, e_vdw_direct,e_elec_direct, atom_voxel );

  //auto step5=std::chrono::system_clock::now();

  //std::chrono::duration<double> d4 = step5 -step4;

  //s: %f seconds\n",d4.count());


  //mprintf("e_elec_self: %f , e_elec_direct: %f, e_vdw6direct: %f \n", e_self, e_direct, e_vdw);

  if (debug_ > 0)
    mprintf("DEBUG: Eself= %20.10f   Erecip= %20.10f   Edirect= %20.10f  Evdw= %20.10f\n",
            e_self, e_recip, e_direct, e_vdw);
  e_vdw += (e_vdw_lr_correction + e_vdw6self + e_vdw6recip);
  t_total_.Stop();
  e_elec = e_self + e_recip + e_direct;
  return 0;
}


/** Electrostatic self energy. This is the cancelling Gaussian plus the "neutralizing plasma". */
double GIST_PME::Self_GIST(double volume, Darray& atom_self) {
  t_self_.Start();
  double d0 = -ew_coeff_ * INVSQRTPI_;
  double ene = SumQ2() * d0;
//  mprintf("DEBUG: d0= %20.10f   ene= %20.10f\n", d0, ene);
  double factor = Constants::PI / (ew_coeff_ * ew_coeff_ * volume);
  double ee_plasma = -0.5 * factor * SumQ() * SumQ();

  for( unsigned int i=0; i< Charge_.size();i++)
  {
    atom_self[i]=Charge_[i]*Charge_[i]*d0 + ee_plasma/Charge_.size(); // distribute the "neutrilizing plasma" to atoms equally
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
double GIST_PME::Recip_ParticleMesh_GIST(Box const& boxIn, MatType& potential)
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
  double erecip = pme_object_.computeERec(0, chargesD, coordsD);
  pme_object_.computePRec(0,chargesD,coordsD,coordsD,1,potential);




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
double GIST_PME::Vdw_Correction_GIST(double volume, Darray& e_vdw_lr_cor) {
  double prefac = Constants::TWOPI / (3.0*volume*cutoff_*cutoff_*cutoff_);
  //mprintf("VDW correction prefac: %.15f \n", prefac);
  double e_vdwr = -prefac * Vdw_Recip_Term();

  //mprintf("Cparam size: %i \n",Cparam_.size());

  //mprintf("volume of the unit cell: %f", volume);


  for ( unsigned int i = 0; i != Cparam_.size(); i++)
  {
    double term(0);

    int v_type=vdw_type_[i];

    //v_type is the vdw_type of atom i, each atom has a atom type

    // atype_vdw_recip_terms_[vdw_type[i]] is the total vdw_recep_term for this v_type

    //N_vdw_type_[v_type] is the total atom number belongs to this v_type



    term = atype_vdw_recip_terms_[v_type] / N_vdw_type_[v_type];

    //mprintf("for i = %i,vdw_type = %i, Number of atoms in this vdw_type_= %i, Total vdw_recip_terms_ for this type= %f,  vdw_recip_term for atom i= %f \n",i,vdw_type_[i],N_vdw_type_[vdw_type_[i]],atype_vdw_recip_terms_[vdw_type_[i]],term);

    e_vdw_lr_cor[i]= -prefac * term ;
    //mprintf("atom e_vdw_lr_cor: %f \n", -prefac* atom_vdw_recip_terms_[i]);
  }

  if (debug_ > 0) mprintf("DEBUG: Vdw correction %20.10f\n", e_vdwr);
  return e_vdwr;
}

/** Direct space routine for GIST. */
double GIST_PME::Direct_GIST(PairList const& PL, double& evdw_out, Darray& e_vdw_direct, Darray& e_elec_direct,
                             Iarray const& atom_voxel)
{
  if (lw_coeff_ > 0.0)
    return Direct_VDW_LJPME_GIST(PL, evdw_out, e_vdw_direct, e_elec_direct);
  else
    return Direct_VDW_LongRangeCorrection_GIST(PL, evdw_out, atom_voxel);
}

/** Direct space calculation with long range VDW correction for GIST. */
double GIST_PME::Direct_VDW_LongRangeCorrection_GIST(PairList const& PL, double& evdw_out,
                                                  Iarray const& atom_voxel)
{
  // e_vdw_direct only count the interaction water molecule involved( sw, ww)
  // e_vdw_all_direct, count all interactions( sw, ww, ss)

  t_direct_.Start();
  double Eelec = 0.0;
  double e_adjust = 0.0;
  double Evdw = 0.0;
  int cidx;

  //mprintf("Entering Direct_VDW_LongrangeCorrection_GIST function \n");
  double* e_vdw_direct  = &(E_vdw_direct_[0][0]);
  double* e_elec_direct = &(E_elec_direct_[0][0]);
  // Pair list loop
# ifdef _OPENMP
  int mythread;
# pragma omp parallel private(cidx, mythread, e_vdw_direct, e_elec_direct) reduction(+: Eelec, Evdw, e_adjust)
  {
  mythread = omp_get_thread_num();
  e_vdw_direct  = &(E_vdw_direct_[mythread][0]);
  e_elec_direct = &(E_elec_direct_[mythread][0]);
# pragma omp for
# endif
//# include "PairListLoop.h"
  //double Evdw_temp(0);   //Eljpme_correction_temp(0), Eljpme_correction_excl_temp(0);
  //double Eelec_temp(0),E_adjust_temp(0);
  
  
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
//#             include "EnergyKernel_Nonbond.h"
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
                 
                e_elec_direct[it0->Idx()] += 0.5 * e_elec;
                e_elec_direct[it1->Idx()] += 0.5 * e_elec;

                int nbindex = NB().GetLJindex(TypeIdx(it0->Idx()),
                                              TypeIdx(it1->Idx()));
                if (nbindex > -1) {
                  double vswitch = SwitchFxn(rij2, cut2_0_, cut2_);
                  NonbondType const& LJ = NB().NBarray()[ nbindex ];
                  double r2    = 1.0 / rij2;
                  double r6    = r2 * r2 * r2;
                  double r12   = r6 * r6;
                  double f12   = LJ.A() * r12;  // A/r^12
                  double f6    = LJ.B() * r6;   // B/r^6
                  double e_vdw = f12 - f6;      // (A/r^12)-(B/r^6)
                  Evdw += (e_vdw * vswitch);

                  // add the vdw energy to e_vdw_direct, divide the energy evenly to each atom

                  e_vdw_direct[it0->Idx()] += 0.5 * e_vdw * vswitch;
                  e_vdw_direct[it1->Idx()] += 0.5 * e_vdw * vswitch;

                  //mprintf("PVDW %8i%8i%20.6f%20.6f\n", ta0+1, ta1+1, e_vdw, r2);
#                 ifdef CPPTRAJ_EKERNEL_LJPME
                  // LJ PME direct space correction
                  double kr2 = lw_coeff_ * lw_coeff_ * rij2;
                  double kr4 = kr2 * kr2;
                  //double kr6 = kr2 * kr4;
                  double expterm = exp(-kr2);
                  double Cij = Cparam_[it0->Idx()] * Cparam_[it1->Idx()];
                  Eljpme_correction += (1.0 - (1.0 +  kr2 + kr4/2.0)*expterm) * r6 * vswitch * Cij;
#                 endif
                }
            }
          } else {
//#           include "EnergyKernel_Adjust.h"

              double adjust = AdjustFxn(q0,q1,sqrt(rij2));

              e_adjust += adjust;

              //add the e_elec to E_elec_direct[it0->Idx] and E_elec_direct[it1-> Idx], 0.5 to avoid double counting

              e_elec_direct[it0->Idx()] += 0.5 * adjust;
              e_elec_direct[it1->Idx()] += 0.5 * adjust;

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
//#               include "EnergyKernel_Nonbond.h"
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

                e_elec_direct[it0->Idx()] += 0.5 * e_elec;
                e_elec_direct[it1->Idx()] += 0.5 * e_elec;

                int nbindex = NB().GetLJindex(TypeIdx(it0->Idx()),
                                              TypeIdx(it1->Idx()));
                if (nbindex > -1) {
                  double vswitch = SwitchFxn(rij2, cut2_0_, cut2_);
                  NonbondType const& LJ = NB().NBarray()[ nbindex ];
                  double r2    = 1.0 / rij2;
                  double r6    = r2 * r2 * r2;
                  double r12   = r6 * r6;
                  double f12   = LJ.A() * r12;  // A/r^12
                  double f6    = LJ.B() * r6;   // B/r^6
                  double e_vdw = f12 - f6;      // (A/r^12)-(B/r^6)
                  Evdw += (e_vdw * vswitch);

                  // add the vdw energy to e_vdw_direct

                  e_vdw_direct[it0->Idx()] += 0.5 * e_vdw * vswitch;
                  e_vdw_direct[it1->Idx()] += 0.5 * e_vdw * vswitch;


                  //mprintf("PVDW %8i%8i%20.6f%20.6f\n", ta0+1, ta1+1, e_vdw, r2);
#                 ifdef CPPTRAJ_EKERNEL_LJPME
                  // LJ PME direct space correction
                  double kr2 = lw_coeff_ * lw_coeff_ * rij2;
                  double kr4 = kr2 * kr2;
                  //double kr6 = kr2 * kr4;
                  double expterm = exp(-kr2);
                  double Cij = Cparam_[it0->Idx()] * Cparam_[it1->Idx()];
                  Eljpme_correction += (1.0 - (1.0 +  kr2 + kr4/2.0)*expterm) * r6 * vswitch * Cij;
#                 endif
                }

              }
            } else {
//#             include "EnergyKernel_Adjust.h"
              double adjust = AdjustFxn(q0,q1,sqrt(rij2));

              e_adjust += adjust;

              e_elec_direct[it0->Idx()] += 0.5 * adjust;
              e_elec_direct[it1->Idx()] += 0.5 * adjust;

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
          } // END loop over neighbor cell atoms
        
        } // END Loop over neighbor cells
        
      } // Loop over thisCell atoms
      


    } // END if thisCell is not empty
  } // Loop over cells

# ifdef _OPENMP
  } // END pragma omp parallel
# endif
  t_direct_.Stop();
# ifdef DEBUG_PAIRLIST
  mprintf("DEBUG: Elec                             = %16.8f\n", Eelec);
  mprintf("DEBUG: Eadjust                          = %16.8f\n", e_adjust);
  mprintf("DEBUG: LJ vdw                           = %16.8f\n", Evdw);
# endif
  evdw_out = Evdw;
  return Eelec + e_adjust;
}

/** Direct space calculation with LJ PME for GIST. */ // TODO enable
double GIST_PME::Direct_VDW_LJPME_GIST(PairList const& PL, double& evdw_out, Darray& e_vdw_direct, Darray& e_elec_direct)
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
