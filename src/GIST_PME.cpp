#ifdef LIBPME
#include "GIST_PME.h"
#include "Frame.h"
#include "helpme_standalone.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
GIST_PME::GIST_PME() {}

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
  //t_total_.Start();
  double volume = frameIn.BoxCrd().CellVolume();

  //auto step0 = std::chrono::system_clock::now();


  double e_self = Self_GIST( volume, e_elec_self ); // decomposed E_elec_self for GIST

  //auto step1 = std::chrono::system_clock::now();

  //std::chrono::duration<double> d1 = step1 -step0;

  //mprintf("Eelec_self takes: %f seconds\n", d1.count());


  double e_vdw_lr_correction;

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
  helpme::Matrix<double> e_potentialD(atom_num,4);

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
  double e_vdw6self, e_vdw6recip;
  if (lw_coeff_ > 0.0) {
    helpme::Matrix<double> vdw_potentialD(atom_num,4);

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
double GIST_PME::Self_GIST(double volume, std::vector<double>& atom_self) {
  //t_self_.Start();
  double d0 = -ew_coeff_ * INVSQRTPI_;
  double ene = sumq2_ * d0;
//  mprintf("DEBUG: d0= %20.10f   ene= %20.10f\n", d0, ene);
  double factor = Constants::PI / (ew_coeff_ * ew_coeff_ * volume);
  double ee_plasma = -0.5 * factor * sumq_ * sumq_;

  for( unsigned int i=0; i< Charge_.size();i++)
  {
    atom_self[i]=Charge_[i]*Charge_[i]*d0 + ee_plasma/Charge_.size(); // distribute the "neutrilizing plasma" to atoms equally
  }

  ene += ee_plasma;
  //t_self_.Stop();
  return ene;
}


#endif /* LIBPME */
