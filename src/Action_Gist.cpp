#include <cmath>
using namespace std;
#include "Action_Gist.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" 
#include "Constants.h" // GASK_KCAL and SMALL

// CONSTRUCTOR
Action_Gist::Action_Gist() :
  CurrentParm_(0),
//  watermodel_(false),
//  useTIP3P_(false),
//  useTIP4P_(false),
//  useTIP4PEW_(false),
//  useTIP5P_(false),
//  useTIP3PFW_(false),
//  useSPCE_(false),
//  useSPCFW_(false),
  doOrder_(false),
  doEij_(false)
{
  BULK_DENS_ = 0;
  gridcntr_[0] = -1;
  gridcntr_[1] = -1;
  gridcntr_[2] = -1;
    
  gridorig_[0] = -1;
  gridorig_[1] = -1;
  gridorig_[2] = -1;
  
  gridspacn_ = 0;
 } 


void Action_Gist::Help() const {
//  mprintf("<watermodel>[{tip3p|tip4p|tip4pew}] [doorder] [doeij] [gridcntr <xval> <yval> <zval>] [griddim <xval> <yval> <zval>] [gridspacn <spaceval>] [out <filename>] \n");
  mprintf("\t[doorder] [doeij] [skipE] [refdens <rdval>] [Temp <tval>] [gridcntr <xval> <yval> <zval>]\n"
          "\t[griddim <xval> <yval> <zval>] [gridspacn <spaceval>]\n"
          "\t[out <filename>]\n");
/*  mprintf("\tGIST needs the specification of the water model being used. Supported water models are: \n");
  mprintf("\ta) TIP3P specified as tip3p. \n");
  mprintf("\tb) TIP4P specified as tip4p. \n");
  mprintf("\tc) TIP4PEW specified as tip4pew. \n");
  mprintf("\td) TIP5P specified as tip5p. \n");
  mprintf("\te) TIP3PFW specified as tip3pfw. \n");
  mprintf("\tf) SPCE specified as spce. \n");
  mprintf("\tg) SPCFW specified as spcfw. \n");
  mprintf("  Calculate GIST between water molecules in selected region \n");
*/}

// Action_Gist::Init()
Action::RetType Action_Gist::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  if (init.TrajComm().Size() > 1) {
    mprinterr("Error: 'gist' action does not work with > 1 thread (%i threads currently).\n",
              init.TrajComm().Size());
    return Action::ERR;
  }
# endif
  if (init.DSL().EnsembleNum() > -1) {
    mprinterr("Error: GIST currently cannot be used in ensemble mode.\n");
    return Action::ERR;
  }
  gist_init_.Start();
  // Get keywords
  // Dataset to store gist results
  datafile_ = actionArgs.GetStringKey("out");
  // Generate the data set name, and hold onto the master data set list
  /*string ds_name = actionArgs.GetStringKey("name");
  ds_name = myDSL_.GenerateDefaultName("GIST");
  // We have 4?? data sets Add them here
  // Now add all of the data sets
  for (int i = 0; i < 4; i++) {
    myDSL_.AddSetAspect(DataSet::DOUBLE, ds_name,
      integerToString(i+1).c_str());
      }
  //  myDSL_.AddSet(DataSet::DOUBLE, ds_name, NULL);
  */  
/*  useTIP3P_ = actionArgs.hasKey("tip3p");
  useTIP4P_ = actionArgs.hasKey("tip4p");
  useTIP4PEW_ = actionArgs.hasKey("tip4pew");
  useTIP5P_ = actionArgs.hasKey("tip5p");
  useTIP3PFW_ = actionArgs.hasKey("tip3pfw");
  useSPCE_ = actionArgs.hasKey("spce");
  useSPCFW_ = actionArgs.hasKey("spcfw");
  if (!useTIP3P_ && !useTIP4P_ && !useTIP4PEW_ && !useTIP5P_ && !useTIP3PFW_ && !useSPCE_ && !useSPCFW_) {
    mprinterr("Init Error: Only water models supported are TIP3P, TIP4P, TIP4PEW, TIP5P, TIP3P/FW, SPC/E, SPC/FW\n");
    return Action::ERR;
  }
*/
  doOrder_ = actionArgs.hasKey("doorder");
  doEij_ = actionArgs.hasKey("doeij");
  skipE_ = actionArgs.hasKey("skipE");
  gridspacn_ = actionArgs.getKeyDouble("gridspacn", 0.50);
  // Set Bulk Energy based on water model
/*  if (useTIP3P_) BULK_E_ = -19.0653;
  if (useTIP4PEW_) BULK_E_ = -22.071;
  if (useTIP4P_) BULK_E_ = -19.71152;
  if (useTIP5P_) BULK_E_ = -19.19174;
  if (useTIP3PFW_) BULK_E_ = -22.7374;
  if (useSPCE_) BULK_E_ = -22.2574;
  if (useSPCFW_) BULK_E_ = -23.7458;
//  if (usePOL3_) BULK_E_ = -22.071;
  mprintf("\tGIST bulk energy: %10.5f\n", BULK_E_);
*/  
  // Set Bulk Density 55.5M
  BULK_DENS_ = actionArgs.getKeyDouble("refdens", -1);
  Temp = actionArgs.getKeyDouble("temp", -1);
  // Grid center
  if ( actionArgs.hasKey("gridcntr") ) {
    gridcntr_[0] = actionArgs.getNextDouble(-1);
    gridcntr_[1] = actionArgs.getNextDouble(-1);
    gridcntr_[2] = actionArgs.getNextDouble(-1);
  } else {
    mprintf("Warning: No grid center values specified, using default\n");
    gridcntr_[0] = 0.0;
    gridcntr_[1] = 0.0;
    gridcntr_[2] = 0.0;
  }
  // Grid dimensions
  griddim_.clear();
  griddim_.resize( 3 );
  if ( actionArgs.hasKey("griddim") ) {
    griddim_[0] = actionArgs.getNextInteger(-1);
    griddim_[1] = actionArgs.getNextInteger(-1);
    griddim_[2] = actionArgs.getNextInteger(-1);
  } else {
    mprintf("Warning: No grid dimension values specified, using default\n");
    griddim_[0] = 40;
    griddim_[1] = 40;
    griddim_[2] = 40;
  }

  InitImaging(true); // Always image

  mprintf("    GIST:\n");
  if(doOrder_)
    mprintf("\tDo Order calculation\n");
  else
    mprintf("\tSkip Order calculation\n");
  if(doEij_)
    mprintf("\tCompute and print water-water Eij matrix\n");
  else
    mprintf("\tSkip water-water Eij matrix\n");
  if (BULK_DENS_ < 0) {
    BULK_DENS_ = 0.0334;
  //BULK_DENS_ = 0.033422885325;
    mprintf("\tNo water reference density specified, using default: %6.4f, equivalent to 1g/cc\n",
            BULK_DENS_);
  } else {
    mprintf("\tWater reference density: %6.4f\n", BULK_DENS_);
    if ( BULK_DENS_ > (0.0334*1.2) )
      mprintf("Warning: water reference density is high, consider using 0.0334 for 1g/cc water density\n");
    else if ( BULK_DENS_ < (0.0334*0.8) )
      mprintf("Warning: water reference density is low, consider using 0.0334 for 1g/cc water density\n");
  }
  if (Temp < 0) {
    Temp = 300.0;
    mprintf("\tNo simulation temperature specified, using default: %6.4f\n", Temp);
  }
  mprintf("\tGIST grid center: %5.3f %5.3f %5.3f\n", gridcntr_[0],gridcntr_[1],gridcntr_[2]);
  mprintf("\tGIST grid dimension: %d %d %d \n", griddim_[0],griddim_[1],griddim_[2]);
  mprintf("\tGIST grid spacing: %5.3f A^3\n", gridspacn_);
  mprintf("\t#Please cite these papers if you use GIST results in a publication:\n"
	  "\t#    Steven Ramsey, Crystal Nguyen, Romelia Salomon-Ferrer, Ross C. Walker, Michael K. Gilson, and Tom Kurtzman J. Comp. Chem. 37 (21) 2016\n"
          "\t#    Crystal Nguyen, Michael K. Gilson, and Tom Young, arXiv:1108.4876v1 (2011)\n"
          "\t#    Crystal N. Nguyen, Tom Kurtzman Young, and Michael K. Gilson,\n"
          "\t#      J. Chem. Phys. 137, 044101 (2012)\n"
          "\t#    Lazaridis, J. Phys. Chem. B 102, 3531–3541 (1998)\n");
  gist_init_.Stop();
  return Action::OK;
}

// Action_Gist::Setup()
/** Set Gist up for this parmtop. Get masks etc.
  */
Action::RetType Action_Gist::Setup(ActionSetup& setup) {
  gist_setup_.Start();
  CurrentParm_ = setup.TopAddress();
  NFRAME_ = 0;
  max_nwat_ = 0;

  MAX_GRID_PT_ = griddim_[0] * griddim_[1] * griddim_[2];
  Vvox_ = gridspacn_*gridspacn_*gridspacn_;
  G_max_x_ = griddim_[0] * gridspacn_ + 1.5 ;
  G_max_y_ = griddim_[1] * gridspacn_ + 1.5 ;
  G_max_z_ = griddim_[2] * gridspacn_ + 1.5 ;
  
  //mprintf("\tGIST Setup: %d %d %d %d %f \n", griddim_[0], griddim_[1], 
  //        griddim_[2], MAX_GRID_PT_, Vvox_);
  mprintf("\tGIST number of voxels: %d, voxel volume: %f A^3\n",  MAX_GRID_PT_, Vvox_);

  // Set up grid origin
  gridorig_[0] = gridcntr_[0] - 0.5*griddim_[0]*gridspacn_;
  gridorig_[1] = gridcntr_[1] - 0.5*griddim_[1]*gridspacn_;
  gridorig_[2] = gridcntr_[2] - 0.5*griddim_[2]*gridspacn_;
  mprintf("\tGIST grid origin: %5.3f %5.3f %5.3f\n", 
          gridorig_[0], gridorig_[1], gridorig_[2]);

  // Set up cumulative energy arrays
  /*  x_.clear();
  x_.resize(5, 0.0);
  y_.clear();
  y_.resize(5, 0.0);
  z_.clear();
  z_.resize(5, 0.0);*/
  wh_evdw_.clear();
  wh_evdw_.resize(MAX_GRID_PT_, 0.0);
  wh_eelec_.clear();
  wh_eelec_.resize(MAX_GRID_PT_, 0.0);
  ww_evdw_.clear();
  ww_evdw_.resize(MAX_GRID_PT_, 0.0);
  ww_eelec_.clear();
  ww_eelec_.resize(MAX_GRID_PT_, 0.0);

  //voxel coords
  grid_x_.clear();    
  grid_x_.resize(MAX_GRID_PT_, 0.0); 
  grid_y_.clear();          
  grid_y_.resize(MAX_GRID_PT_, 0.0);
  grid_z_.clear();          
  grid_z_.resize(MAX_GRID_PT_, 0.0); 


  // get the actual voxel coordinates
  voxel_ = 0;
  for (int i = 0; i < griddim_[0]; ++i) {
    for (int j = 0; j < griddim_[1]; ++j) {
      for (int k = 0; k < griddim_[2]; ++k) {
        grid_x_[voxel_] = Xcrd(i);
        grid_y_[voxel_] = Ycrd(j);
        grid_z_[voxel_] = Zcrd(k);
        voxel_++;
      }
    }
  }

  Esw_dens_.clear();
  Esw_dens_.resize(MAX_GRID_PT_, 0.0);
  Esw_norm_.clear();
  Esw_norm_.resize(MAX_GRID_PT_, 0.0);
  Eww_norm_.clear();
  Eww_norm_.resize(MAX_GRID_PT_, 0.0);
  Eww_dens_.clear();
  Eww_dens_.resize(MAX_GRID_PT_, 0.0);

  if(doEij_) {
    ww_Eij_.clear();
    ww_Eij_.resize(MAX_GRID_PT_);
    for(int i = 1; i < MAX_GRID_PT_; i++) ww_Eij_[i].resize(i);
    
    //CN: need to initialize ww_Eij_ to 0.0 but not Euler angles
    for (int a=1; a<MAX_GRID_PT_; a++)
      for (int l=0; l<a; l++) ww_Eij_[a][l]=0.0;  
  }
  x_vox_.clear();
  x_vox_.resize(MAX_GRID_PT_);
  y_vox_.clear();
  y_vox_.resize(MAX_GRID_PT_);
  z_vox_.clear();
  z_vox_.resize(MAX_GRID_PT_);
  q0_vox_.clear();
  q0_vox_.resize(MAX_GRID_PT_);
  q1_vox_.clear();
  q1_vox_.resize(MAX_GRID_PT_);
  q2_vox_.clear();
  q2_vox_.resize(MAX_GRID_PT_);
  q3_vox_.clear();
  q3_vox_.resize(MAX_GRID_PT_);
  
  
  dTStrans_dens_.clear();
  dTStrans_dens_.resize(MAX_GRID_PT_, 0.0);
  dTStrans_norm_.clear();
  dTStrans_norm_.resize(MAX_GRID_PT_, 0.0); 
  dTSorient_dens_.clear();
  dTSorient_dens_.resize(MAX_GRID_PT_, 0.0);
  dTSorient_norm_.clear();
  dTSorient_norm_.resize(MAX_GRID_PT_, 0.0);
  dTSsix_norm_.clear();
  dTSsix_norm_.resize(MAX_GRID_PT_, 0.0);
  dTSsix_dens_.clear();
  dTSsix_dens_.resize(MAX_GRID_PT_, 0.0);

  nwat_.clear();
  nwat_.resize(MAX_GRID_PT_, 0);
  nH_.clear();
  nH_.resize(MAX_GRID_PT_, 0);
  nw_angle_.clear();
  nw_angle_.resize(MAX_GRID_PT_, 0);
  dens_.clear();
  dens_.resize(MAX_GRID_PT_, 0.0);
  g_.clear();
  g_.resize(MAX_GRID_PT_, 0.0);
  gH_.clear();
  gH_.resize(MAX_GRID_PT_, 0.0);
  dipolex_.clear();
  dipolex_.resize(MAX_GRID_PT_, 0.0);
  dipoley_.clear();
  dipoley_.resize(MAX_GRID_PT_, 0.0);
  dipolez_.clear();
  dipolez_.resize(MAX_GRID_PT_, 0.0);
  neighbor_.clear();
  neighbor_.resize(MAX_GRID_PT_, 0.0);
  neighbor_dens_.clear();
  neighbor_dens_.resize(MAX_GRID_PT_, 0.0);
  neighbor_norm_.clear();
  neighbor_norm_.resize(MAX_GRID_PT_, 0.0);
  qtet_.clear();
  qtet_.resize(MAX_GRID_PT_, 0.0);
  pol_.clear();
  pol_.resize(MAX_GRID_PT_, 0.0);

  gridwat_.clear();
  gridwat_.resize( setup.Top().Nsolvent() );

  // We need box info
  if (setup.CoordInfo().TrajBox().Type() == Box::NOBOX) {
    mprinterr("Error: Gist: Must have explicit solvent with periodic boundaries!");
    return Action::ERR;
  }
  SetupImaging( setup.CoordInfo().TrajBox().Type() );

  resnum_ = 0;
  voxel_ = 0;
  gist_setup_.Stop();
  return Action::OK;  
}

// Action_Gist::DoAction()
Action::RetType Action_Gist::DoAction(int frameNum, ActionFrame& frm) {
  NFRAME_ ++;
//  if (NFRAME_==1) mprintf("GIST Action \n");

  // Simulation box length - assign here because it can vary for npt simulation
  //Lx_ = frm.Frm().BoxCrd().BoxX();
  //Ly_ = frm.Frm().BoxCrd().BoxY();
  //Lz_ = frm.Frm().BoxCrd().BoxZ();
//  if (NFRAME_==1) mprintf("GIST Action box length: %f %f %f \n", Lx_, Ly_, Lz_);
  
  int solventMolecules = CurrentParm_->Nsolvent();
  resnum_ = 0;
  voxel_ = 0;
  resindex1_ = 0;
  for (solvmol_ = CurrentParm_->MolStart();
       solvmol_ != CurrentParm_->MolEnd(); ++solvmol_)
  {
    if (skipE_) {
      resindex1_++;
      if (!solvmol_->IsSolvent()) continue;
      gist_grid_.Start();
      Grid( frm.Frm() );
      gist_grid_.Stop();
      voxel_ = gridwat_[resnum_];
      resnum_++;
      if (voxel_ >= MAX_GRID_PT_) continue;
      gist_euler_.Start();
      EulerAngle( frm.Frm() );
      gist_euler_.Stop();
      gist_dipole_.Start();
      Dipole( frm.Frm() );
      gist_dipole_.Stop();
    }
    else {
      resindex1_++;
      if (!solvmol_->IsSolvent()) continue;
      gist_grid_.Start();
      Grid( frm.Frm() );
      gist_grid_.Stop();
      voxel_ = gridwat_[resnum_];
      resnum_++;
      gist_nonbond_.Start();
      NonbondEnergy( frm.Frm() );
      gist_nonbond_.Stop();
      if (voxel_ >= MAX_GRID_PT_) continue;
      gist_euler_.Start();
      EulerAngle( frm.Frm() );
      gist_euler_.Stop();
      gist_dipole_.Start();
      Dipole( frm.Frm() );
      gist_dipole_.Stop();
    }
  }
  if(doOrder_) Order( frm.Frm() );
  
  //Debug
//  if (NFRAME_==1) mprintf("GIST  DoAction:  Found %d solvent residues \n", resnum_);
  if (solventMolecules != resnum_) {
    mprinterr("GIST  DoAction Error: Number of solvent molecules don't match %d %d\n", solventMolecules, resnum_);
  }
  
  return Action::OK;
}

// Action_Gist::NonbondEnergy()
void Action_Gist::NonbondEnergy(Frame const& frameIn) {
  double rij2, rij, r2, r6, r12, f12, f6, e_vdw, e_elec;
  int satom, satom2, atom1, atom2;
  
  int  voxel2 = 0;
  double q1, q2;
  
  // Setup imaging info
  Matrix_3x3 ucell, recip;
  if (ImagingEnabled())
    frameIn.BoxCrd().ToRecip(ucell, recip);

  // Inner loop has both solute and solvent
  resnum2_=0;
  resindex2_ = 1;
  // skip if water2 has index larger than water1 so that every pair is only evaluated once
  solvmol2_ = CurrentParm_->MolStart();
  for (resindex2_=1; resindex2_<resindex1_; resindex2_++)
  {    
    if (!solvmol2_->IsSolvent()) {
      // Outer loop is not water, break inner loop if water 1 is outside the grid
      if (voxel_ >= MAX_GRID_PT_) {
        ++solvmol2_;
        continue;
      }
    } else { 
      // Inner loop is water
      voxel2 = gridwat_[resnum2_];
      resnum2_++;
      // skip if both waters are outside the grid
      if (voxel_>=MAX_GRID_PT_ && voxel2>=MAX_GRID_PT_) {
        ++solvmol2_;
        continue;
      }
    }
      
    // Loop over all solvent atoms of water 1
    atom1=0;
    for (satom = solvmol_->BeginAtom(); satom < solvmol_->EndAtom(); ++satom)
    {
      // Set up coord index for this atom
      const double* XYZ =  frameIn.XYZ( satom );
      atom2=0;
      for (satom2 = solvmol2_->BeginAtom(); satom2 < solvmol2_->EndAtom(); ++satom2)
      {    
        // Set up coord index for this atom
        const double* XYZ2 = frameIn.XYZ( satom2 );
        // Calculate the vector pointing from atom2 to atom1
        rij2 = DIST2(XYZ, XYZ2, ImageType(), frameIn.BoxCrd(), ucell, recip);
        rij = sqrt(rij2);
        // LJ energy
        NonbondType const& LJ = CurrentParm_->GetLJparam(satom, satom2); 
        r2    = 1 / rij2;
        r6    = r2 * r2 * r2;
        r12   = r6 * r6;
        f12   = LJ.A() * r12;  // A/r^12
        f6    = LJ.B() * r6;   // B/r^6
        e_vdw = f12 - f6;     // (A/r^12)-(B/r^6)
        // LJ Force 
        // Coulomb energy 
        q1 = (*CurrentParm_)[satom].Charge() * Constants::ELECTOAMBER;
        q2 = (*CurrentParm_)[satom2].Charge() * Constants::ELECTOAMBER;
        e_elec = (q1*q2/rij);
        if (!solvmol2_->IsSolvent()) {
          // solute-solvent interaction
          wh_evdw_[voxel_] +=  e_vdw;
          wh_eelec_[voxel_] += e_elec;
        } else {
          // solvent-solvent interaction, need to compute for all waters,
          // even those outside the grid but only one water needs to be 
          // inside the grid. 
          //mprintf("DEBUG0: v1= %i v2= %i EVV %i %i Vdw= %f Elec= %f\n", voxel_, voxel2, satom, satom2, e_vdw, e_elec);
          if (voxel_<MAX_GRID_PT_) {
            ww_evdw_[voxel_] +=  e_vdw;
            ww_eelec_[voxel_] += e_elec;
            // Store the water neighbor using only O-O distance
            if (atom2==0 && atom1==0 && rij<3.5)
              neighbor_[voxel_] += 1.0;
          }
          // CN: only store Eij[voxel1][voxel2] if both voxels lie on the grid.
          if (voxel2<MAX_GRID_PT_) {
            ww_evdw_[voxel2] +=  e_vdw;
            ww_eelec_[voxel2] += e_elec;
            // Store the water neighbor using only O-O distance
            if (atom2==0 && atom1==0 && rij<3.5)
              neighbor_[voxel2] += 1.0;
            if (doEij_ && (voxel_<MAX_GRID_PT_)) {
              if (voxel_>voxel2) {
                ww_Eij_[voxel_][voxel2] += e_vdw;
                ww_Eij_[voxel_][voxel2] += e_elec;
              } else {
                ww_Eij_[voxel2][voxel_] += e_vdw;
                ww_Eij_[voxel2][voxel_] += e_elec;
              }
            } //print Eij && voxel<MAX_GRID_PT_
          }
        }//IF is solvent
        atom2++;
      } // END Inner loop ALL atoms
      atom1++;
    } // END Outer loop solvent atoms
    ++solvmol2_;
  }  // END Inner loop ALL molecules
  //if (voxel_ < MAX_GRID_PT_) mprintf("DEBUG0: atom %i voxel %i VV evdw=%f eelec=%f\n", solvmol_->BeginAtom(), voxel_, ww_evdw_[voxel_], ww_eelec_[voxel_]);
}

// Action_Gist::Grid()
void Action_Gist::Grid(Frame const& frameIn) {
  int  i, gridindex[3], nH;
  Vec3 comp,  atom_coord;
  i = solvmol_->BeginAtom();

  gridwat_[resnum_] = MAX_GRID_PT_ + 1;
  atom_coord = Vec3(frameIn.XYZ(i));
  // get the components of the water vector
  comp = Vec3(atom_coord) - Vec3(gridorig_);
  nH=0;
  //If Oxygen is far from grid, 1.5A or more in any durection, skip calculation
  if (comp[0] <= G_max_x_ && comp[1] <= G_max_y_ && comp[2] <= G_max_z_ && 
      comp[0] >= -1.5     && comp[1] >= -1.5     && comp[2] >= -1.5 )
  {
    //if (comp[0]<= G_max_x || comp[1]<= G_max_y || comp[2]<= G_max_z ||
    //    comp[0]>= -1.5 || comp[1]>= -1.5 || comp[2]>= -1.5 ) {
    //Water is at most 1.5A away from grid, so we need to check for H even if O is outside grid
    nH=2;
    
    //O is inside grid only if comp is >=0
    if (comp[0]>=0 && comp[1]>=0 && comp[2]>=0 ){
      comp /= gridspacn_;
      gridindex[0] = (int) comp[0];
      gridindex[1] = (int) comp[1];
      gridindex[2] = (int) comp[2];
      
      if ((gridindex[0]<griddim_[0]) && (gridindex[1]<griddim_[1]) && (gridindex[2]<griddim_[2]))
      {
        // this water belongs to grid point gridindex[0], gridindex[1], gridindex[2]
        voxel_ = (gridindex[0]*griddim_[1] + gridindex[1])*griddim_[2] + gridindex[2];
        gridwat_[resnum_] = voxel_;
        //mprintf("DEBUG0: Water atom %i voxel %i\n", i, voxel_);
        nwat_[voxel_]++;
        if (max_nwat_ < nwat_[voxel_]) max_nwat_ = nwat_[voxel_];
      }
    }
    
    // evaluate hydrogen atoms
    for (int a=1; a<=nH; a++) {
      atom_coord = Vec3(frameIn.XYZ(i+a));
      comp = Vec3(atom_coord) - Vec3(gridorig_);
      if (comp[0]<0 || comp[1]<0 || comp[2]<0) continue;
      comp /= gridspacn_;
      gridindex[0] = (int) comp[0];
      gridindex[1] = (int) comp[1];
      gridindex[2] = (int) comp[2];
      if ((gridindex[0]<griddim_[0]) && 
          (gridindex[1]<griddim_[1]) && 
          (gridindex[2]<griddim_[2]))
      {
        voxel_ = (gridindex[0]*griddim_[1] + gridindex[1])*griddim_[2] + gridindex[2];
        nH_[voxel_]++;
      }
    } 
  }
}

// Action_Gist::EulerAngle()
void Action_Gist::EulerAngle(Frame const& frameIn) {
  //if (NFRAME_==1) mprintf("GIST Euler Angles \n");
  Vec3 x_lab, y_lab, z_lab, O_wat, H1_wat, H2_wat, x_wat, y_wat, z_wat, node, v, ar1, ar2, ar3;
  //double dp; double w1, w2, w3, w4, w5, x1, x2, x3, x4, x5, y1, y2, y3, y4, y5, z1, z2, z3, z4, z5;
  double w1, w2, w3, w4, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
  //double dp1, dp2, dp3, sign; Vec3 sar, H_temp, H_temp2;
  double dp1, dp2, sign; Vec3 sar, H_temp, H_temp2;

  int i = solvmol_->BeginAtom();
  O_wat = Vec3(frameIn.XYZ(i));
  x_vox_[voxel_].push_back(O_wat[0]); y_vox_[voxel_].push_back(O_wat[1]); z_vox_[voxel_].push_back(O_wat[2]);
  H1_wat = Vec3(frameIn.XYZ(i+1)) - O_wat;
  H2_wat = Vec3(frameIn.XYZ(i+2)) - O_wat;
  //Steve
  O_wat = O_wat - O_wat;

  // make sure the first three atoms are oxygen followed by two hydrogen
  if ((*CurrentParm_)[i].Element() != Atom::OXYGEN) {
    mprintf("Warning: GIST: First coordinates do not belong to oxygen atom (%s)\n",
            (*CurrentParm_)[i].ElementName());
  }
  if ((*CurrentParm_)[i+1].Element() != Atom::HYDROGEN ||
      (*CurrentParm_)[i+2].Element() != Atom::HYDROGEN)
  {
    mprintf("Warning: GIST: second and third coordinates do not belong to hydrogen atoms (%s, %s)\n",
            (*CurrentParm_)[i+1].ElementName(), (*CurrentParm_)[i+2].ElementName());
  }

  // Define lab frame of reference
  x_lab[0]=1.0; x_lab[1]=0;   x_lab[2]=0;
  y_lab[0]=0;   y_lab[1]=1.0; y_lab[2]=0;
  z_lab[0]=0;   z_lab[1]=0;   z_lab[2]=1.0;

  //Steve
  H1_wat.Normalize();
  H2_wat.Normalize();
  ar1 = H1_wat.Cross(x_lab);
  ar1.Normalize();
  dp1 = x_lab*H1_wat;
  theta_ = acos(dp1);
  sar = H1_wat.Cross(x_lab);
  sign = sar*H1_wat;
  if (sign > 0) theta_/=2;
  else theta_/=-2;
  w1 = cos(theta_);
  x1 = ar1[0]*sin(theta_);
  y1 = ar1[1]*sin(theta_);
  z1 = ar1[2]*sin(theta_);
  w2 = w1; x2 = x1; y2 = y1; z2 = z1;
  H_temp[0] = ((w2*w2+x2*x2)-(y2*y2+z2*z2))*H1_wat[0];
  H_temp[0] = (2*(x2*y2 - w2*z2)*H1_wat[1]) + H_temp[0];
  H_temp[0] = (2*(x2*z2-w2*y2)*H1_wat[2]) + H_temp[0];

  H_temp[1] = 2*(x2*y2 - w2*z2)* H1_wat[0];
  H_temp[1] = ((w2*w2-x2*x2+y2*y2-z2*z2)*H1_wat[1]) + H_temp[1];
  H_temp[1] = (2*(y2*z2+w2*x2)*H1_wat[2]) +H_temp[1];

  H_temp[2] = 2*(x2*z2+w2*y2) *H1_wat[0];
  H_temp[2] = (2*(y2*z2-w2*x2)*H1_wat[1]) + H_temp[2];
  H_temp[2] = ((w2*w2-x2*x2-y2*y2+z2*z2)*H1_wat[2]) + H_temp[2];

  H1_wat = H_temp;

  H_temp2[0] = ((w2*w2+x2*x2)-(y2*y2+z2*z2))*H2_wat[0];
  H_temp2[0] = (2*(x2*y2 + w2*z2)*H2_wat[1]) + H_temp2[0];
  H_temp2[0] = (2*(x2*z2-w2*y2)+H2_wat[2]) +H_temp2[0];

  H_temp2[1] = 2*(x2*y2 - w2*z2) *H2_wat[0];
  H_temp2[1] = ((w2*w2-x2*x2+y2*y2-z2*z2)*H2_wat[1]) +H_temp2[1];
  H_temp2[1] = (2*(y2*z2+w2*x2)*H2_wat[2]) +H_temp2[1];

  H_temp2[2] = 2*(x2*z2+w2*y2)*H2_wat[0];
  H_temp2[2] = (2*(y2*z2-w2*x2)*H2_wat[1]) +H_temp2[2];
  H_temp2[2] = ((w2*w2-x2*x2-y2*y2+z2*z2)*H2_wat[2]) + H_temp2[2];

  H2_wat = H_temp2;

  ar2 = H_temp.Cross(H_temp2);
  ar2.Normalize();
  dp2 = ar2*z_lab;
  theta_ = acos(dp2);

  sar = ar2.Cross(z_lab);
  sign = sar*H_temp;

  if (sign < 0) theta_/=2;
  else theta_/=-2;

  w3 = cos(theta_);
  x3 = x_lab[0]*sin(theta_);
  y3 = x_lab[1]*sin(theta_);
  z3 = x_lab[2]*sin(theta_);

  w4 = w1*w3 - x1*x3 - y1*y3 - z1*z3;
  x4 = w1*x3 + x1*w3 + y1*z3 - z1*y3;
  y4 = w1*y3 - x1*z3 + y1*w3 + z1*x3;
  z4 = w1*z3 + x1*y3 - y1*x3 + z1*w3;

  q0_vox_[voxel_].push_back(w4);
  q1_vox_[voxel_].push_back(x4);
  q2_vox_[voxel_].push_back(y4);
  q3_vox_[voxel_].push_back(z4);
  //mprintf("DEBUG0: wxyz4= %g %g %g %g\n", w4, x4, y4, z4);
  nw_angle_[voxel_]++;
  //}
  //else mprintf("%i: gimbal lock problem, two z_wat paralell\n", resnum-1);
} 

// Action_Gist::Dipole()
void Action_Gist::Dipole(Frame const& frameIn) {
  
  //if (NFRAME_==1) mprintf("GIST Dipole \n");
  double dipolar_vector[3], charge;
  int satom;

  dipolar_vector[0] = 0.0;
  dipolar_vector[1] = 0.0;
  dipolar_vector[2] = 0.0;
  // Loop over solvent atoms
  for (satom = solvmol_->BeginAtom(); satom < solvmol_->EndAtom(); ++satom)
  {
    const double* XYZ = frameIn.XYZ( satom );
    // Calculate dipole vector. The oxygen of the solvent is used to 
    // assign the voxel index to the water.
    // NOTE: the total charge on the solvent should be neutral for this 
    //       to have any meaning
    charge = (*CurrentParm_)[satom].Charge();
    //      mprintf("%i %f %f %f %f\n", resnum-1, charge, XYZ[0], XYZ[1], XYZ[2]);
    dipolar_vector[0] += (charge * XYZ[0]);
    dipolar_vector[1] += (charge * XYZ[1]);
    dipolar_vector[2] += (charge * XYZ[2]);
  }
  //mprintf("DEBUG0: voxel %i dipole %f %f %f\n", voxel_, dipolar_vector[0], dipolar_vector[1], dipolar_vector[2]);
  dipolex_[voxel_] += dipolar_vector[0];
  dipoley_[voxel_] += dipolar_vector[1];
  dipolez_[voxel_] += dipolar_vector[2];
}

// Action_Gist::Order() 
void Action_Gist::Order(Frame const& frameIn) {
//  if (NFRAME_==1) mprintf("GIST Order Parameter \n");
  int i;
  double cos, sum, r1, r2, r3, r4, rij2, x[5], y[5], z[5];
  Vec3 neighbor1(0.0), neighbor2(0.0), neighbor3(0.0), neighbor4(0.0);
  Vec3 O_wat1, O_wat2, O_wat3, v1, v2;
  resnum_=0;

  for (solvmol_ = CurrentParm_->MolStart();
       solvmol_ != CurrentParm_->MolEnd(); ++solvmol_)
  {
    if (!solvmol_->IsSolvent()) continue;

    // obtain 4 closest neighbors for every water
    resnum_++;
    voxel_ = gridwat_[resnum_-1];
    if (voxel_>=MAX_GRID_PT_) continue;
    // assume that oxygen is the first atom
    i = solvmol_->BeginAtom();
    O_wat1 = Vec3(frameIn.XYZ( i ));

    r1=1000; r2=1000; r3=1000; r4=1000; resnum2_=0;
    // Can't make into triangular matrix
    for (solvmol2_ = CurrentParm_->MolStart();
         solvmol2_ != CurrentParm_->MolEnd(); ++solvmol2_)
    {
      if (!solvmol2_->IsSolvent()) continue;
      resnum2_++;
      if (resnum_ == resnum2_) continue;
      i = solvmol2_->BeginAtom();
      O_wat2 = Vec3(frameIn.XYZ( i ));      
      rij2 = DIST2_NoImage(O_wat1, O_wat2);
      if (rij2<r1) {
        r4 = r3;
        r3 = r2;
        r2 = r1;
        r1 = rij2;
        neighbor4 = neighbor3;
        neighbor3 = neighbor2;
        neighbor2 = neighbor1;
        neighbor1 = O_wat2;
      }
      else if (rij2<r2) {
        r4 = r3;
        r3 = r2;
        r2 = rij2;
        neighbor4 = neighbor3;
        neighbor3 = neighbor2;
        neighbor2 = O_wat2;
      }
      else if (rij2<r3) {
        r4 = r3;
        r3 = rij2;
        neighbor4 = neighbor3;
        neighbor3 = O_wat2;       
      }
      else if (rij2<r4) {
        r4 = rij2;
        neighbor4 = O_wat2;
     }        
    }
    x[1]=neighbor1[0]; y[1]=neighbor1[1]; z[1]=neighbor1[2];
    x[2]=neighbor2[0]; y[2]=neighbor2[1]; z[2]=neighbor2[2];
    x[3]=neighbor3[0]; y[3]=neighbor3[1]; z[3]=neighbor3[2];
    x[4]=neighbor4[0]; y[4]=neighbor4[1]; z[4]=neighbor4[2];    
    // Compute the tetrahedral order parameter
    sum=0;
    for (int mol1=1; mol1<=3; mol1++) {
      for (int mol2=mol1+1; mol2<=4; mol2++) {
        O_wat2[0] = x[mol1];
        O_wat2[1] = y[mol1];
        O_wat2[2] = z[mol1];
        O_wat3[0] = x[mol2];
        O_wat3[1] = y[mol2];
        O_wat3[2] = z[mol2];
        v1 = O_wat2 - O_wat1;
        v2 = O_wat3 - O_wat1;    
        r1 = v1.Magnitude2();
        r2 = v2.Magnitude2();
        cos = (v1*( v2))/sqrt(r1*r2);
        sum += (cos + 1.0/3)*(cos + 1.0/3);
      }
    }
    qtet_[voxel_] += (1.0 - (3.0/8)*sum);
/*    double dbl = (1.0 - (3.0/8)*sum)
 *    if (dbl<-3.0 || dbl>1.0) {
      mprintf("BAD! voxel=%d, q=%9.5f\n", voxel_, dbl);
    }*/
  }
}

void Action_Gist::Print() {
  gist_print_.Start();
  int addx = griddim_[2]*griddim_[1];
  int addy = griddim_[2];
  int addz = 1;
  int subx = -1*griddim_[2]*griddim_[1];
  int suby = -1*griddim_[2];
  int subz = -1;
  // Implement NN to compute orientational entropy for each voxel
  //double NNr, rx, ry, rz, rR, dbl;
  double NNr, rR, dbl;
  double dTSs = 0.0, dTSst = 0.0;
  int nwtt = 0;
  float NNs = 10000;
  float ds = 0;
  float dx = 0, dy = 0, dz = 0, dd = 0, NNd = 10000;
  double dTSo = 0, dTSot = 0, dTSt = 0, dTStt = 0;
  int nwts = 0;

  double dTSorienttot_ = 0; //NFRAME_ /= 8;
  for (int gr_pt = 0; gr_pt < MAX_GRID_PT_; gr_pt++) {
    //int numplane = gr_pt/(griddim_[1]*griddim_[2]); int nwj = 0;
    dTSorient_dens_[gr_pt]=0;
    dTSorient_norm_[gr_pt]=0;
    if (nw_angle_[gr_pt] != nwat_[gr_pt])
      mprintf("DEBUG: voxel %i nw_angle_=%i nwat_=%i\n", gr_pt, nw_angle_[gr_pt], nwat_[gr_pt]);
    int nwtot = nw_angle_[gr_pt];
    //mprintf("DEBUG0: %i nw_total %i\n", gr_pt, nwtot);
    int bound = 0;
    nwtt += nwtot;
    if (nwtot <= 1) continue;
    for (int n = 0; n < nwtot; n++) {
      NNr = 10000;
      NNs = 10000;
      ds = 0;
      NNd = 10000;
      dd = 0;
      for (int l = 0; l < nwtot; l++) {
        if (l == n) continue;
        rR = 2 * acos(  q0_vox_[gr_pt][l] * q0_vox_[gr_pt][n]
                      + q1_vox_[gr_pt][l] * q1_vox_[gr_pt][n]
                      + q2_vox_[gr_pt][l] * q2_vox_[gr_pt][n]
                      + q3_vox_[gr_pt][l] * q3_vox_[gr_pt][n] );
        //mprintf("DEBUG0: %g\n", rR);
        if (rR>0 && rR < NNr) NNr = rR;
      }

      if (bound == 1) {
        dbl = 0;
        dTSorient_norm_[gr_pt] += dbl;
        continue;
      } else if (NNr < 9999 && NNr > 0 && NNs > 0) {
        dbl = log(NNr*NNr*NNr*nwtot / (3.0*Constants::TWOPI));
        //mprintf("DEBUG0: dbl %f\n", dbl);
        dTSorient_norm_[gr_pt] += dbl;
        dTSo += dbl;
      }
    } // END loop over nwtot
    //NFRAME_ *= 0.5;
    //mprintf("DEBUG0: dTSorient_norm %f\n", dTSorient_norm_[gr_pt]);
    dTSorient_norm_[gr_pt] = Constants::GASK_KCAL * Temp * // FIXME what is this constant?
                             ((dTSorient_norm_[gr_pt]/nwtot)+0.5772156649);
    dTSorient_dens_[gr_pt] = dTSorient_norm_[gr_pt] * nwat_[gr_pt] / (NFRAME_ * Vvox_);
    dTSorienttot_ += dTSorient_dens_[gr_pt];
    //mprintf("DEBUG0: %f\n", dTSorienttot_);
  } // END loop over MAX_GRID_PT_
  dTSorienttot_ *= Vvox_;
  mprintf("Maximum number of waters found in one voxel for %d frames = %d\n", NFRAME_, max_nwat_);
  mprintf("Total referenced orientational entropy of the grid: dTSorient = %9.5f kcal/mol, Nf=%d\n",
          dTSorienttot_, NFRAME_);

  // Compute translational entropy for each voxel
  double dTStranstot_ = 0.0;
    for (int a = 0; a < MAX_GRID_PT_; a++) {
      int numplane = a / (griddim_[1]*griddim_[2]);
      int nwj = 0;
      dens_[a] = 1.0*nwat_[a]/(NFRAME_*Vvox_);
      g_[a] = dens_[a] / BULK_DENS_;
      gH_[a] = 1.0 * nH_[a] / (NFRAME_*Vvox_*2*BULK_DENS_);
      int nwi = 0, bound = 0;

      //first do own voxel
      nwi = nwat_[a];
      for (int i = 0; i < nwi; i++) {
        NNd = 10000;
        bound = 0;
        NNs = 10000;
        ds = 0;
        for (int j = 0; j < nwi; j++) {
          if ( j!= i) {
            //cout << "doing self: " << a << endl;
            dx = x_vox_[a][i] - x_vox_[a][j];
            dy = y_vox_[a][i] - y_vox_[a][j];
            dz = z_vox_[a][i] - z_vox_[a][j];
            dd = dx*dx+dy*dy+dz*dz;
            if (dd < NNd && dd > 0) { NNd = dd; }
            rR = 2 * acos( q0_vox_[a][i]*q0_vox_[a][j]
                          +q1_vox_[a][i]*q1_vox_[a][j]
                          +q2_vox_[a][i]*q2_vox_[a][j]
                          +q3_vox_[a][i]*q3_vox_[a][j] );
            ds = rR*rR + dd;
            if (ds < NNs && ds > 0) { NNs = ds; }
          }
        } // END loop over nwi
        //mprintf("DEBUG0: self NNd=%f NNs=%f\n", NNd, NNs);

        //if (a+addz > MAX_GRID_PT_ || a+addz < 0) {throw exc;}
        if (griddim_[2] == 0 || (a%griddim_[2] == griddim_[2]-1))
          bound = 1;
        else {
          //cout << "doing addz: " << a << endl;
          //mprintf("else bound, addz: %d\n", a);
          nwj = nwat_[a+addz];
          for (int j = 0; j < nwj; j++) {
            dx = x_vox_[a][i] - x_vox_[a+addz][j];
            dy = y_vox_[a][i] - y_vox_[a+addz][j];
            dz = z_vox_[a][i] - z_vox_[a+addz][j];
            dd = dx*dx+dy*dy+dz*dz;
            if (dd < NNd && dd > 0) { NNd = dd; }
            rR = 2 * acos( q0_vox_[a][i]*q0_vox_[a+addz][j]
                          +q1_vox_[a][i]*q1_vox_[a+addz][j]
                          +q2_vox_[a][i]*q2_vox_[a+addz][j]
                          +q3_vox_[a][i]*q3_vox_[a+addz][j]);
            ds = rR*rR + dd;
            if (ds < NNs && ds > 0) { NNs = ds; }
          }
        }

        //if (a+addy > MAX_GRID_PT_ || a+addy < 0) {throw exc;}
        if ((griddim_[2] == 0 || griddim_[1]-1 == 0) ||
            (a%(griddim_[2]*(griddim_[1]-1)+(numplane*griddim_[2]*griddim_[1])) < griddim_[2]))
          bound = 1;
        else {
          //cout << "doing addy: " << a << endl;
          //mprintf("else bound, addy: %d\n", a);
          nwj = nwat_[a+addy];
          for (int j = 0; j < nwj; j++) {
            dx = x_vox_[a][i] - x_vox_[a+addy][j];
            dy = y_vox_[a][i] - y_vox_[a+addy][j];
            dz = z_vox_[a][i] - z_vox_[a+addy][j];
            dd = dx*dx+dy*dy+dz*dz;
            if (dd < NNd && dd > 0) { NNd = dd; }
            rR = 2 * acos( q0_vox_[a][i]*q0_vox_[a+addy][j]
                          +q1_vox_[a][i]*q1_vox_[a+addy][j]
                          +q2_vox_[a][i]*q2_vox_[a+addy][j]
                          +q3_vox_[a][i]*q3_vox_[a+addy][j]);
            ds = rR*rR + dd;
            if (ds < NNs && ds > 0) {NNs = ds;}
          }
        }

        //if (a+addx > MAX_GRID_PT_ || a+addx < 0) {throw exc;}
        if (a >= griddim_[2]*griddim_[1] * (griddim_[0]-1) &&
            a <  griddim_[2]*griddim_[1] *  griddim_[0]      )
          bound = 1;
        else {
          //cout << "doing addx: " << a << endl;
          //mprintf("else bound, addx: %d\n", a);
          nwj = nwat_[a+addx];
          for (int j = 0; j < nwj; j++) {
            dx = x_vox_[a][i] - x_vox_[a+addx][j];
            dy = y_vox_[a][i] - y_vox_[a+addx][j];
            dz = z_vox_[a][i] - z_vox_[a+addx][j];
            dd = dx*dx+dy*dy+dz*dz;
            if (dd < NNd && dd > 0) { NNd = dd; }
            rR = 2 * acos( q0_vox_[a][i]*q0_vox_[a+addx][j]
                          +q1_vox_[a][i]*q1_vox_[a+addx][j]
                          +q2_vox_[a][i]*q2_vox_[a+addx][j]
                          +q3_vox_[a][i]*q3_vox_[a+addx][j]);
            ds = rR*rR + dd;
            if (ds < NNs && ds > 0) { NNs = ds; }
          }
        }

        //if (a+subz > MAX_GRID_PT_ || a+subz < 0) {throw exc;}
        if (griddim_[2] == 0 || a%griddim_[2] == 0)
          bound = 1;
        else {
          //cout << "doing subz: " << a << endl;
          //mprintf("else bound, subz: %d\n", a);
          nwj = nwat_[a+subz];
          for (int j = 0; j < nwj; j++) {
            dx = x_vox_[a][i] - x_vox_[a+subz][j];
            dy = y_vox_[a][i] - y_vox_[a+subz][j];
            dz = z_vox_[a][i] - z_vox_[a+subz][j];
            dd = dx*dx+dy*dy+dz*dz;
            if (dd < NNd && dd > 0) {NNd = dd;}
            rR = 2 * acos( q0_vox_[a][i]*q0_vox_[a+subz][j]
                          +q1_vox_[a][i]*q1_vox_[a+subz][j]
                          +q2_vox_[a][i]*q2_vox_[a+subz][j]
                          +q3_vox_[a][i]*q3_vox_[a+subz][j]);
            ds = rR*rR + dd;
            if (ds < NNs && ds > 0) {NNs = ds;}
          }
        }

        //if (a+suby > MAX_GRID_PT_ || a+suby < 0) {throw exc;}
        if ( (griddim_[2] == 0 || griddim_[1] == 0) ||
             (a%(griddim_[2]*griddim_[1]) < griddim_[2]) )
          bound = 1;
        else {
          //cout << "doing suby: " << a << endl;
          //mprintf("else bound, suby: %d\n", a);
          nwj = nwat_[a+suby];
          for (int j = 0; j < nwj; j++) {
            dx = x_vox_[a][i] - x_vox_[a+suby][j];
            dy = y_vox_[a][i] - y_vox_[a+suby][j];
            dz = z_vox_[a][i] - z_vox_[a+suby][j];
            dd = dx*dx+dy*dy+dz*dz;
            if (dd < NNd && dd > 0) {NNd = dd;}
            rR = 2 * acos( q0_vox_[a][i]*q0_vox_[a+suby][j]
                          +q1_vox_[a][i]*q1_vox_[a+suby][j]
                          +q2_vox_[a][i]*q2_vox_[a+suby][j]
                          +q3_vox_[a][i]*q3_vox_[a+suby][j]);
            ds = rR*rR + dd;
            if (ds < NNs && ds > 0) {NNs = ds;}
          }
        }

        //if (a+subx > MAX_GRID_PT_ || a+subx < 0) {throw exc;}
        if ((griddim_[2] == 0 || griddim_[1] == 0) ||
            (a >= 0 && a < griddim_[2]*griddim_[1]))
          bound = 1;
        else {
          //cout << "doing subx: " << a << endl;
          //mprintf("else bound, subx: %d\n", a);
          nwj = nwat_[a+subx];
          for (int j = 0; j < nwj; j++) {
              dx = x_vox_[a][i] - x_vox_[a+subx][j];
              dy = y_vox_[a][i] - y_vox_[a+subx][j];
              dz = z_vox_[a][i] - z_vox_[a+subx][j];
              dd = dx*dx+dy*dy+dz*dz;
              if (dd < NNd && dd > 0) {NNd = dd;}
              rR = 2 * acos( q0_vox_[a][i]*q0_vox_[a+subx][j]
                            +q1_vox_[a][i]*q1_vox_[a+subx][j]
                            +q2_vox_[a][i]*q2_vox_[a+subx][j]
                            +q3_vox_[a][i]*q3_vox_[a+subx][j]);
              ds = rR*rR + dd;
              if (ds < NNs && ds > 0) {NNs = ds;}
          }
        }

        //if (a+addz+addy > MAX_GRID_PT_ || a+addz+addy < 0) {throw exc;}
        if ((griddim_[2] == 0 || griddim_[1]-1 == 0) ||
            (a%griddim_[2] == griddim_[2]-1) ||
            (a%(griddim_[2]*(griddim_[1]-1)+(numplane*griddim_[2]*griddim_[1])) < griddim_[2]))
          bound = 1;
        else {
          //cout << "doing addz addy: " << a << endl;
          //mprintf("else bound, addz+addy: %d\n", a);
          nwj = nwat_[a+addz+addy];
          for (int j = 0; j < nwj; j++) {
            dx = x_vox_[a][i] - x_vox_[a+addz+addy][j];
            dy = y_vox_[a][i] - y_vox_[a+addz+addy][j];
            dz = z_vox_[a][i] - z_vox_[a+addz+addy][j];
            dd = dx*dx+dy*dy+dz*dz;
            if (dd < NNd && dd > 0) {NNd = dd;}
            rR = 2 * acos( q0_vox_[a][i]*q0_vox_[a+addz+addy][j]
                          +q1_vox_[a][i]*q1_vox_[a+addz+addy][j]
                          +q2_vox_[a][i]*q2_vox_[a+addz+addy][j]
                          +q3_vox_[a][i]*q3_vox_[a+addz+addy][j]);
            ds = rR*rR + dd;
            if (ds < NNs && ds > 0) {NNs = ds;}
          }
        }

        //if (a+addz+suby > MAX_GRID_PT_ || a+addz+suby < 0) {throw exc;}
        if ((griddim_[1] == 0 || griddim_[2] == 0) ||
            (a%griddim_[2] == griddim_[2]-1) || 
            (a%(griddim_[2]*griddim_[1]) < griddim_[2]))
          bound = 1;
        else {
          //cout << "doing addz suby: " << a << endl;
          //mprintf("else bound, addz+suby: %d\n", a);
          nwj = nwat_[a+addz+suby];
          for (int j = 0; j < nwj; j++) {
            dx = x_vox_[a][i] - x_vox_[a+addz+suby][j];
            dy = y_vox_[a][i] - y_vox_[a+addz+suby][j];
            dz = z_vox_[a][i] - z_vox_[a+addz+suby][j];
            dd = dx*dx+dy*dy+dz*dz;
            if (dd < NNd && dd > 0) {NNd = dd;}
            rR = 2 * acos( q0_vox_[a][i]*q0_vox_[a+addz+suby][j]
                          +q1_vox_[a][i]*q1_vox_[a+addz+suby][j]
                          +q2_vox_[a][i]*q2_vox_[a+addz+suby][j]
                          +q3_vox_[a][i]*q3_vox_[a+addz+suby][j]);
            ds = rR*rR + dd;
            if (ds < NNs && ds > 0) {NNs = ds;}
          }
        }

        //if (a+subz+addy > MAX_GRID_PT_ || a+subz+addy < 0) {throw exc;}
        if ((griddim_[2] == 0 || griddim_[1]-1 == 0) ||
            (a%griddim_[2] == 0) ||
            (a%(griddim_[2]*(griddim_[1]-1)+(numplane*griddim_[2]*griddim_[1]))< griddim_[2]))
          bound = 1;
        else {
          //cout << "doing subz addy: " << a << endl;
          //mprintf("else bound, subz+addy: %d\n", a);
          nwj = nwat_[a+subz+addy];
          for (int j = 0; j < nwj; j++) {
            dx = x_vox_[a][i] - x_vox_[a+subz+addy][j];
            dy = y_vox_[a][i] - y_vox_[a+subz+addy][j];
            dz = z_vox_[a][i] - z_vox_[a+subz+addy][j];
            dd = dx*dx+dy*dy+dz*dz;
            if (dd < NNd && dd > 0) {NNd = dd;}
            rR = 2 * acos( q0_vox_[a][i]*q0_vox_[a+subz+addy][j]
                          +q1_vox_[a][i]*q1_vox_[a+subz+addy][j]
                          +q2_vox_[a][i]*q2_vox_[a+subz+addy][j]
                          +q3_vox_[a][i]*q3_vox_[a+subz+addy][j]);
            ds = rR*rR + dd;
            if (ds < NNs && ds > 0) {NNs = ds;}
          }
        }

        //if (a+subz+suby > MAX_GRID_PT_ || a+subz+suby < 0) {throw exc;}
        if ((griddim_[2] == 0 || griddim_[1] == 0) ||
            (a%griddim_[2] == 0) ||
            (a%(griddim_[2]*griddim_[1]) < griddim_[2]))
          bound = 1;
        else {
          //cout << "doing subz suby: " << a << endl;
          //mprintf("else bound, subz+suby: %d\n", a);
          nwj = nwat_[a+subz+suby];
          for (int j = 0; j < nwj; j++) {
            dx = x_vox_[a][i] - x_vox_[a+subz+suby][j];
            dy = y_vox_[a][i] - y_vox_[a+subz+suby][j];
            dz = z_vox_[a][i] - z_vox_[a+subz+suby][j];
            dd = dx*dx+dy*dy+dz*dz;
            if (dd < NNd && dd > 0) {NNd = dd;}
            rR = 2 * acos( q0_vox_[a][i]*q0_vox_[a+subz+suby][j]
                          +q1_vox_[a][i]*q1_vox_[a+subz+suby][j]
                          +q2_vox_[a][i]*q2_vox_[a+subz+suby][j]
                          +q3_vox_[a][i]*q3_vox_[a+subz+suby][j]);
            ds = rR*rR + dd; if (ds < NNs && ds > 0) {NNs = ds;}
          }
        }

        //if (a+addz+addx > MAX_GRID_PT_ || a+addz+addx < 0) {throw exc;}
        if ((griddim_[2] == 0) ||
            (a%griddim_[2] == griddim_[2]-1) ||
            (a >= griddim_[2]*griddim_[1]*(griddim_[0]-1) &&
             a < griddim_[2]*griddim_[1]*griddim_[0]))
          bound = 1;
        else {
          //cout << "doing addz addx: " << a << endl;
          //mprintf("else bound, addz+addx: %d\n", a);
          nwj = nwat_[a+addz+addx];
          for (int j = 0; j < nwj; j++) {
            dx = x_vox_[a][i] - x_vox_[a+addz+addx][j];
            dy = y_vox_[a][i] - y_vox_[a+addz+addx][j];
            dz = z_vox_[a][i] - z_vox_[a+addz+addx][j];
            dd = dx*dx+dy*dy+dz*dz;
            if (dd < NNd && dd > 0) {NNd = dd;}
            rR = 2 * acos( q0_vox_[a][i]*q0_vox_[a+addz+addx][j]
                          +q1_vox_[a][i]*q1_vox_[a+addz+addx][j]
                          +q2_vox_[a][i]*q2_vox_[a+addz+addx][j]
                          +q3_vox_[a][i]*q3_vox_[a+addz+addx][j]);
            ds = rR*rR + dd;
            if (ds < NNs && ds > 0) {NNs = ds;}
          }
        }

        //if (a+addz+subx > MAX_GRID_PT_ || a+addz+subx < 0) {throw exc;}
        if ((griddim_[2] == 0) ||
            (a%griddim_[2] == griddim_[2]-1) ||
            (a >= 0 && a < griddim_[2]*griddim_[1]))
          bound = 1;
        else {
          //cout << "doing addz subx: " << a << endl;
          //mprintf("else bound, addz+subx: %d\n", a);
          nwj = nwat_[a+addz+subx];
          for (int j = 0; j < nwj; j++) {
            dx = x_vox_[a][i] - x_vox_[a+addz+subx][j];
            dy = y_vox_[a][i] - y_vox_[a+addz+subx][j];
            dz = z_vox_[a][i] - z_vox_[a+addz+subx][j];
            dd = dx*dx+dy*dy+dz*dz;
            if (dd < NNd && dd > 0) {NNd = dd;}
            rR = 2 * acos( q0_vox_[a][i]*q0_vox_[a+addz+subx][j]
                          +q1_vox_[a][i]*q1_vox_[a+addz+subx][j]
                          +q2_vox_[a][i]*q2_vox_[a+addz+subx][j]
                          +q3_vox_[a][i]*q3_vox_[a+addz+subx][j]);
            ds = rR*rR + dd;
            if (ds < NNs && ds > 0) {NNs = ds;}
          }
        }

        //if (a+subz+addx > MAX_GRID_PT_ || a+subz+addx < 0) {throw exc;}
        if ((griddim_[2] == 0) ||
            (a%griddim_[2] == 0)||
            (a >= griddim_[2]*griddim_[1]*(griddim_[0]-1) &&
             a < griddim_[2]*griddim_[1]*griddim_[0]))
          bound = 1;
        else {
          //cout << "doing subz addx: " << a << endl;
          //mprintf("else bound, subz+addx: %d\n", a);
          nwj = nwat_[a+subz+addx];
          for (int j = 0; j < nwj; j++) {
            dx = x_vox_[a][i] - x_vox_[a+subz+addx][j];
            dy = y_vox_[a][i] - y_vox_[a+subz+addx][j];
            dz = z_vox_[a][i] - z_vox_[a+subz+addx][j];
            dd = dx*dx+dy*dy+dz*dz;
            if (dd < NNd && dd > 0) {NNd = dd;}
            rR = 2 * acos( q0_vox_[a][i]*q0_vox_[a+subz+addx][j]
                          +q1_vox_[a][i]*q1_vox_[a+subz+addx][j]
                          +q2_vox_[a][i]*q2_vox_[a+subz+addx][j]
                          +q3_vox_[a][i]*q3_vox_[a+subz+addx][j]);
            ds = rR*rR + dd;
            if (ds < NNs && ds > 0) {NNs = ds;}
          }
        }

        //if (a+subz+subx > MAX_GRID_PT_ || a+subz+subx < 0) {throw exc;}
        if ((griddim_[2] == 0) ||
            (a%griddim_[2] == 0) ||
            (a >= 0 && a < griddim_[2]*griddim_[1]))
          bound = 1;
        else {
          //cout << "doing subz subx: " << a << endl;
          //mprintf("else bound, subz+subx: %d\n", a);
          nwj = nwat_[a+subz+subx];
          for (int j = 0; j < nwj; j++) {
            dx = x_vox_[a][i] - x_vox_[a+subz+subx][j];
            dy = y_vox_[a][i] - y_vox_[a+subz+subx][j];
            dz = z_vox_[a][i] - z_vox_[a+subz+subx][j];
            dd = dx*dx+dy*dy+dz*dz;
            if (dd < NNd && dd > 0) {NNd = dd;}
            rR = 2 * acos( q0_vox_[a][i]*q0_vox_[a+subz+subx][j]
                          +q1_vox_[a][i]*q1_vox_[a+subz+subx][j]
                          +q2_vox_[a][i]*q2_vox_[a+subz+subx][j]
                          +q3_vox_[a][i]*q3_vox_[a+subz+subx][j]);
            ds = rR*rR + dd;
            if (ds < NNs && ds > 0) {NNs = ds;}
          }
        }

        //if (a+addy+addx > MAX_GRID_PT_ || a+addy+addx < 0) {throw exc;}
        if ((griddim_[2] == 0 || griddim_[1]-1 == 0) ||
            (a%(griddim_[2]*(griddim_[1]-1)+
                (numplane*griddim_[2]*griddim_[1])) < griddim_[2]) ||
            (a >= griddim_[2]*griddim_[1]*(griddim_[0]-1) &&
             a < griddim_[2]*griddim_[1]*griddim_[0]))
          bound = 1;
        else {
          //cout << "doing addy addx: " << a << endl;
          //mprintf("else bound, addy+addx: %d\n", a);
          nwj = nwat_[a+addy+addx];
          for (int j = 0; j < nwj; j++) {
            dx = x_vox_[a][i] - x_vox_[a+addy+addx][j];
            dy = y_vox_[a][i] - y_vox_[a+addy+addx][j];
            dz = z_vox_[a][i] - z_vox_[a+addy+addx][j];
            dd = dx*dx+dy*dy+dz*dz;
            if (dd < NNd && dd > 0) {NNd = dd;}
            rR = 2 * acos( q0_vox_[a][i]*q0_vox_[a+addy+addx][j]
                          +q1_vox_[a][i]*q1_vox_[a+addy+addx][j]
                          +q2_vox_[a][i]*q2_vox_[a+addy+addx][j]
                          +q3_vox_[a][i]*q3_vox_[a+addy+addx][j]);
            ds = rR*rR + dd;
            if (ds < NNs && ds > 0) {NNs = ds;}
          }
        }

        //if (a+addy+subx > MAX_GRID_PT_ || a+addy+subx < 0) {throw exc;}
        if ((griddim_[2] == 0 || griddim_[1]-1 == 0) ||
            (a%(griddim_[2]*(griddim_[1]-1)+(numplane*griddim_[2]*griddim_[1]))< griddim_[2]) ||
            (a >= 0 && a < griddim_[2]*griddim_[1]))
          bound = 1;
        else {
          //cout << "doing addy subx: " << a << endl;
          //mprintf("else bound, addy+subx: %d\n", a);
          nwj = nwat_[a+addy+subx];
          for (int j = 0; j < nwj; j++) {
            dx = x_vox_[a][i] - x_vox_[a+addy+subx][j];
            dy = y_vox_[a][i] - y_vox_[a+addy+subx][j];
            dz = z_vox_[a][i] - z_vox_[a+addy+subx][j];
            dd = dx*dx+dy*dy+dz*dz;
            if (dd < NNd && dd > 0) {NNd = dd;}
            rR = 2 * acos( q0_vox_[a][i]*q0_vox_[a+addy+subx][j]
                          +q1_vox_[a][i]*q1_vox_[a+addy+subx][j]
                          +q2_vox_[a][i]*q2_vox_[a+addy+subx][j]
                          +q3_vox_[a][i]*q3_vox_[a+addy+subx][j]);
            ds = rR*rR + dd;
            if (ds < NNs && ds > 0) {NNs = ds;}
          }
        }

        //if (a+suby+addx > MAX_GRID_PT_ || a+suby+addx < 0) {throw exc;}
        if ((griddim_[2] == 0 || griddim_[1] == 0) ||
            (a%(griddim_[2]*griddim_[1]) < griddim_[2]) ||
            (a >= griddim_[2]*griddim_[1]*(griddim_[0]-1) &&
             a < griddim_[2]*griddim_[1]*griddim_[0]))
          bound = 1;
        else {
          //cout << "doing suby addx: " << a << endl;
          //mprintf("else bound, suby+addx: %d\n", a);
          nwj = nwat_[a+suby+addx];
          for (int j = 0; j < nwj; j++) {
            dx = x_vox_[a][i] - x_vox_[a+suby+addx][j];
            dy = y_vox_[a][i] - y_vox_[a+suby+addx][j];
            dz = z_vox_[a][i] - z_vox_[a+suby+addx][j];
            dd = dx*dx+dy*dy+dz*dz;
            if (dd < NNd && dd > 0) {NNd = dd;}
            rR = 2 * acos( q0_vox_[a][i]*q0_vox_[a+suby+addx][j]
                          +q1_vox_[a][i]*q1_vox_[a+suby+addx][j]
                          +q2_vox_[a][i]*q2_vox_[a+suby+addx][j]
                          +q3_vox_[a][i]*q3_vox_[a+suby+addx][j]);
            ds = rR*rR + dd;
            if (ds < NNs && ds > 0) {NNs = ds;}
          }
        }

        //if (a+suby+subx > MAX_GRID_PT_ || a+suby+subx < 0) {throw exc;}
        if ((griddim_[2] == 0 || griddim_[1] == 0) ||
            (a%(griddim_[2]*griddim_[1]) < griddim_[2]) ||
            (a >= 0 && a < griddim_[2]*griddim_[1]))
          bound = 1;
        else {
          //cout << "doing suby subx: " << a << endl;
          //mprintf("else bound, suby+subx: %d\n", a);
          nwj = nwat_[a+suby+subx];
          for (int j = 0; j < nwj; j++) {
            dx = x_vox_[a][i] - x_vox_[a+suby+subx][j];
            dy = y_vox_[a][i] - y_vox_[a+suby+subx][j];
            dz = z_vox_[a][i] - z_vox_[a+suby+subx][j];
            dd = dx*dx+dy*dy+dz*dz;
            if (dd < NNd && dd > 0) {NNd = dd;}
            rR = 2 * acos( q0_vox_[a][i]*q0_vox_[a+suby+subx][j]
                          +q1_vox_[a][i]*q1_vox_[a+suby+subx][j]
                          +q2_vox_[a][i]*q2_vox_[a+suby+subx][j]
                          +q3_vox_[a][i]*q3_vox_[a+suby+subx][j]);
            ds = rR*rR + dd;
            if (ds < NNs && ds > 0) {NNs = ds;}
          }
        }

        NNd = sqrt(NNd);
        NNs = sqrt(NNs);
        if (bound == 1) {
          dbl = 0;
          dTStrans_norm_[a] += dbl;
          continue;
        }// dTSsix_norm_[a] += dbl; continue;}
        else if (NNd < 3 && NNd > 0/*NNd < 9999 && NNd > 0*/) {
          //cout << "calc dbl: " << a << endl;
          dbl = log((NNd*NNd*NNd*NFRAME_*4*Constants::PI*BULK_DENS_)/3);
          dTStrans_norm_[a] += dbl;
          dTSt += dbl;
          dbl = log((NNs*NNs*NNs*NNs*NNs*NNs*NFRAME_*Constants::PI*BULK_DENS_)/48);
          dTSsix_norm_[a] += dbl;
          dTSs += dbl;
          //mprintf("DEBUG0: dbl=%f NNs=%f\n", dbl, NNs);
        }
      }
      if (dTStrans_norm_[a] != 0) {
        //cout << "doing norm: " << a << endl;
        nwts += nwi;
        dTStrans_norm_[a] = Constants::GASK_KCAL*Temp*((dTStrans_norm_[a]/nwi)+0.5772156649);
        dTSsix_norm_[a] = Constants::GASK_KCAL*Temp*((dTSsix_norm_[a]/nwi)+0.5772156649);
      }
      //cout << "doing dens: " << a << endl;
      dTStrans_dens_[a] = dTStrans_norm_[a]*nwat_[a]/(NFRAME_*Vvox_);
      dTSsix_dens_[a] = dTSsix_norm_[a]*nwat_[a]/(NFRAME_*Vvox_);
      dTStranstot_ += dTStrans_dens_[a];
    }

    dTStranstot_ *= Vvox_;
    dTSst = Constants::GASK_KCAL*Temp*((dTSs/nwts) + 0.5772156649);
    dTSot = Constants::GASK_KCAL*Temp*((dTSo/nwtt) + 0.5772156649);
    dTStt = Constants::GASK_KCAL*Temp*((dTSt/nwts) + 0.5772156649);
    mprintf("watcount in vol = %d\n", nwtt);
    mprintf("watcount in subvol = %d\n", nwts);
    mprintf("Total referenced translational entropy of the grid:"
            " dTStrans = %9.5f kcal/mol, Nf=%d\n",
            dTStranstot_, NFRAME_);
    mprintf("Total 6d if all one vox: %9.5f kcal/mol\n", dTSst);
    mprintf("Total t if all one vox: %9.5f kcal/mol\n", dTStt);
    mprintf("Total o if all one vox: %9.5f kcal/mol\n", dTSot);

  // Compute average voxel energy
  double Eswtot = 0.0;
  double Ewwtot = 0.0;
  for (int a=0; a<MAX_GRID_PT_; a++) {
    //mprintf("DEBUG0: VV vdw=%f elec=%f\n", ww_evdw_[a],ww_eelec_[a]);
    if (nwat_[a]>1) {
       Esw_dens_[a] = (wh_evdw_[a]+wh_eelec_[a])/(NFRAME_*Vvox_);
       Esw_norm_[a] = (wh_evdw_[a]+wh_eelec_[a])/nwat_[a];
       Eww_dens_[a] = (ww_evdw_[a]+ww_eelec_[a])/(2*NFRAME_*Vvox_);
       Eww_norm_[a] = (ww_evdw_[a]+ww_eelec_[a])/(2*nwat_[a]);
       Eswtot += Esw_dens_[a];
       Ewwtot += Eww_dens_[a];
    } else {
       Esw_dens_[a]=0; Esw_norm_[a]=0; Eww_norm_[a]=0; Eww_dens_[a]=0;
    }
    // Compute the average number of water neighbor, average order parameter, and average dipole density
    if (nwat_[a]>0) {
      qtet_[a] /= nwat_[a];
      neighbor_norm_[a] = 1.0*neighbor_[a]/nwat_[a];
    }
    neighbor_dens_[a] = 1.0*neighbor_[a]/(NFRAME_*Vvox_);
    dipolex_[a] /= (0.20822678*NFRAME_*Vvox_);
    dipoley_[a] /= (0.20822678*NFRAME_*Vvox_);
    dipolez_[a] /= (0.20822678*NFRAME_*Vvox_);
    pol_[a] = sqrt(dipolex_[a]*dipolex_[a] + dipoley_[a]*dipoley_[a] + dipolez_[a]*dipolez_[a]);
  }
  Eswtot *= Vvox_;
  Ewwtot *= Vvox_;
  mprintf("Total water-solute energy of the grid: Esw = %9.5f kcal/mol\n", Eswtot);
  mprintf("Total unreferenced water-water energy of the grid: Eww = %9.5f kcal/mol\n", Ewwtot);

  // Print the gist info file
  // Print the energy data
  /*if (!datafile_.empty()) {
  // Now write the data file with all of the GIST energies
  DataFile dfl;
  ArgList dummy;
  dfl.SetupDatafile(datafile_, dummy, 0);
  for (int i = 0; i < myDSL_.size(); i++) {
  dfl.AddSet(myDSL_[i]);
  }

  dfl.Write();

  }*/
  //stored as float
  PrintDX("gist-gO.dx", g_);
  PrintDX("gist-gH.dx", gH_);
  PrintDX("gist-Esw-dens.dx", Esw_dens_);
  PrintDX("gist-Eww-dens.dx", Eww_dens_);
  PrintDX("gist-dTStrans-dens.dx", dTStrans_dens_);
  //PrintDX("gist-dTStrans-norm.dx", dTStrans_norm_);
  PrintDX("gist-dTSorient-dens.dx", dTSorient_dens_);
  PrintDX("gist-dTSsix-dens.dx", dTSsix_dens_);
  PrintDX("gist-neighbor-norm.dx", neighbor_norm_);
  PrintDX("gist-dipole-dens.dx", pol_);
  //stored as doubles
  PrintDX_double("gist-order-norm.dx", qtet_);
  PrintDX_double("gist-dipolex-dens.dx", dipolex_);
  PrintDX_double("gist-dipoley-dens.dx", dipoley_);
  PrintDX_double("gist-dipolez-dens.dx", dipolez_);

  if (!datafile_.empty())
    PrintOutput(datafile_);
  else
    PrintOutput("gist-output.dat");
  gist_print_.Stop();
  double total = gist_grid_.Total() + gist_nonbond_.Total() +
                 gist_euler_.Total() + gist_dipole_.Total() +
                 gist_init_.Total() + gist_setup_.Total() +
                 gist_print_.Total();
  mprintf("\tGIST timings:\n");
  gist_init_.WriteTiming(1,    "Init: ", total);
  gist_setup_.WriteTiming(1,   "Setup:", total);
  gist_grid_.WriteTiming(2,    "Grid:   ", total);
  gist_nonbond_.WriteTiming(2, "Nonbond:", total);
  gist_euler_.WriteTiming(2,   "Euler:  ", total);
  gist_dipole_.WriteTiming(2,  "Dipole: ", total);
  gist_print_.WriteTiming(1,   "Print:", total);
  mprintf("TIME:\tTotal: %.4f s\n", total);
}

// Print GIST data in dx format
void Action_Gist::PrintDX_double(string const& filename, std::vector<double>& data)
{
  CpptrajFile outfile;
  if (outfile.OpenWrite(filename)) {
    mprinterr("Print Error: Could not open OpenDX output file.\n");
    return;
  }
  // Print the OpenDX header
  outfile.Printf("object 1 class gridpositions counts %d %d %d\n",
                 griddim_[0], griddim_[1], griddim_[2]);
  outfile.Printf("origin %lg %lg %lg\n", gridorig_[0], gridorig_[1], gridorig_[2]);
  outfile.Printf("delta %lg 0 0\n", gridspacn_);
  outfile.Printf("delta 0 %lg 0\n", gridspacn_);
  outfile.Printf("delta 0 0 %lg\n", gridspacn_);
  outfile.Printf("object 2 class gridconnections counts %d %d %d\n",
                 griddim_[0], griddim_[1], griddim_[2]);
  outfile.Printf(
    "object 3 class array type double rank 0 items %d data follows\n",
    MAX_GRID_PT_);

  // Now print out the data. It is already in row-major form (z-axis changes
  // fastest), so no need to do any kind of data adjustment
  for (int i = 0; i < MAX_GRID_PT_ - 2; i += 3)
    outfile.Printf("%g %g %g\n", data[i], data[i+1], data[i+2]);
  // Print out any points we may have missed
  switch (MAX_GRID_PT_ % 3) {
    case 2: outfile.Printf("%g %g\n", data[MAX_GRID_PT_-2], data[MAX_GRID_PT_-1]); break;
    case 1: outfile.Printf("%g\n", data[MAX_GRID_PT_-1]); break;
  }

  outfile.CloseFile();
}

// Print GIST data in dx format
void Action_Gist::PrintDX(string const& filename, std::vector<float>& data)
{
  CpptrajFile outfile;
  if (outfile.OpenWrite(filename)) {
    mprinterr("Print Error: Could not open OpenDX output file.\n");
    return;
  }
  // Print the OpenDX header
  outfile.Printf("object 1 class gridpositions counts %d %d %d\n",
                 griddim_[0], griddim_[1], griddim_[2]);
  outfile.Printf("origin %lg %lg %lg\n", gridorig_[0], gridorig_[1], gridorig_[2]);
  outfile.Printf("delta %lg 0 0\n", gridspacn_);
  outfile.Printf("delta 0 %lg 0\n", gridspacn_);
  outfile.Printf("delta 0 0 %lg\n", gridspacn_);
  outfile.Printf("object 2 class gridconnections counts %d %d %d\n",
                 griddim_[0], griddim_[1], griddim_[2]);
  outfile.Printf(
    "object 3 class array type float rank 0 items %d data follows\n",
    MAX_GRID_PT_);

  // Now print out the data. It is already in row-major form (z-axis changes
  // fastest), so no need to do any kind of data adjustment
  for (int i = 0; i < MAX_GRID_PT_ - 2; i += 3)
    outfile.Printf("%g %g %g\n", data[i], data[i+1], data[i+2]);
  // Print out any points we may have missed
  switch (MAX_GRID_PT_ % 3) {
    case 2: outfile.Printf("%g %g\n", data[MAX_GRID_PT_-2], data[MAX_GRID_PT_-1]); break;
    case 1: outfile.Printf("%g\n", data[MAX_GRID_PT_-1]); break;
  }

  outfile.CloseFile();
}

// Print GIST data in dx format
void Action_Gist::PrintOutput(string const& filename)
{
  CpptrajFile outfile;
  if (outfile.OpenWrite(filename)) {
    mprinterr("Print Error: Could not open GISToutput file.\n");
    return;
  }
  //ofstream myfile;
  //  myfile.open("gist-output.dat");  
  //myfile.open(filename);  
  outfile.Printf("GIST Output, information printed per voxel\n");
  outfile.Printf("voxel xcoord ycoord zcoord population g_O g_H ");
  outfile.Printf("dTStrans-dens(kcal/mol/A^3) dTStrans-norm(kcal/mol) dTSorient-dens(kcal/mol/A^3) dTSorient-norm(kcal/mol) dTSsix-dens(kcal/mol/A^3) dTSsix-norm (kcal/mol) ");
  outfile.Printf("Esw-dens(kcal/mol/A^3) Esw-norm(kcal/mol) ");
  outfile.Printf("Eww-dens(kcal/mol/A^3) Eww-norm-unref(kcal/mol) ");
  outfile.Printf("Dipole_x-dens(D/A^3) Dipole_y-dens(D/A^3) Dipole_z-dens(D/A^3) Dipole-dens(D/A^3) neighbor-dens(1/A^3) neighbor-norm order-norm\n");
  // Now print out the data. 
  for (int i=0; i<MAX_GRID_PT_; i++){
    outfile.Printf( "%d %g %g %g %d %g %g ",i , grid_x_[i] , grid_y_[i], grid_z_[i], nwat_[i] , g_[i], gH_[i] );
    outfile.Printf( "%g %g %g %g %g %g ",dTStrans_dens_[i], dTStrans_norm_[i], dTSorient_dens_[i] , dTSorient_norm_[i] , dTSsix_dens_[i] , dTSsix_norm_[i] );
    outfile.Printf( "%g %g ",Esw_dens_[i], Esw_norm_[i] );
    outfile.Printf( "%g %g ",Eww_dens_[i] , Eww_norm_[i] );
    outfile.Printf( "%g %g %g %g ",dipolex_[i] , dipoley_[i] , dipolez_[i] , pol_[i] );
    outfile.Printf( "%g %g %g \n",neighbor_dens_[i] , neighbor_norm_[i] , qtet_[i]);
    }
  outfile.CloseFile();

  //double Eijtot=0;
  if(doEij_) {
    if (outfile.OpenWrite("Eww_ij.dat") == 0) {
      double dbl;
      for (int a=1; a < MAX_GRID_PT_; a++) {
        for (int l=0; l<a; l++) {
          dbl = ww_Eij_[a][l];
          if (dbl != 0) {
            dbl /= (NFRAME_*2);
           outfile.Printf("%10d %10d %12.5E\n", a, l, dbl);
           //Eijtot += dbl;
           }
        }
      }
      outfile.CloseFile();
    } else
      mprinterr("Error: Could not open 'Eww_ij.dat' for writing.\n");
  }
  //Eijtot *= 2.0;
  // Debug Eww_ij
  //mprintf("\tTotal grid energy: %9.5f\n", Eijtot);
}
