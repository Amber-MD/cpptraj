// Gist 
#include <cmath>
using namespace std;
#include "Action_Gist.h"
#include "CpptrajFile.h"
#include "CpptrajStdio.h"
#include "DataSet_integer.h"
#include "Box.h"
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


void Action_Gist::Help() {
//  mprintf("<watermodel>[{tip3p|tip4p|tip4pew}] [doorder] [doeij] [gridcntr <xval> <yval> <zval>] [griddim <xval> <yval> <zval>] [gridspacn <spaceval>] [out <filename>] \n");
  mprintf("\t[doorder] [doeij] [refdens <rdval>] [gridcntr <xval> <yval> <zval>]\n"
          "\t[griddim <xval> <yval> <zval>] [gridspacn <spaceval>] [out <filename>]\n");
/*  mprintf("\tGIST needs the specification of the water model being used. Supported water models are: \n");
  mprintf("\ta) TIP3P specified as tip3p. \n");
  mprintf("\tb) TIP4P specified as tip4p. \n");
  mprintf("\tc) TIP4PEW specified as tip4pew. \n");
  mprintf("\td) TIP5P specified as tip5p. \n");
  mprintf("\te) TIP3PFW specified as tip3pfw. \n");
  mprintf("\tf) SPCE specified as spce. \n");
  mprintf("\tg) SPCFW specified as spcfw. \n");
*/}

/*static double diffclock(clock_t clock1,clock_t clock2)
{
  double diffticks=clock1-clock2;
  double diffms=(diffticks*1000)/CLOCKS_PER_SEC;
  return diffms;
}*/


// Action_Gist::Init()
Action::RetType Action_Gist::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                                  DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  if (DSL->EnsembleNum() > -1) {
    mprinterr("Error: GIST currently cannot be used in ensemble mode.\n");
    return Action::ERR;
  }
  mprintf("    GIST:\n");
  // Get keywords
  
//  gist_t_begin=clock();
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
  if(doOrder_){
    mprintf("\tGIST do Order calculation \n");
  }
  else{
    mprintf("\tGIST skip Order calculation \n");
  }

  doEij_ = actionArgs.hasKey("doeij");
  if(doEij_){
    mprintf("\tGIST print water-water Eij matrix \n");
  }
  else{
    mprintf("\tGIST skip water-water Eij matrix \n");
  }

  // Set Bulk Energy based on water model
/*  if (useTIP3P_) BULK_E_ = -19.0653;
  if (useTIP4PEW_) BULK_E_ = -22.071;
  if (useTIP4P_) BULK_E_ = -19.71152;
  if (useTIP5P_) BULK_E_ = -19.19174;
  if (useTIP3PFW_) BULK_E_ = -22.7374;
  if (useSPCE_) BULK_E_ = -22.2574;
  if (useSPCFW_) BULK_E_ = -23.7458;
//  if (usePOL3_) BULK_E_ = -22.071;
  mprintf("\tGIST bulk energy: %10.5f\n", BULK_E_);*/

  // Set Bulk Density 55.5M if no input found
  if ( actionArgs.hasKey("refdens") ){
    BULK_DENS_ = actionArgs.getNextDouble(-1);
    mprintf("\tGIST water reference density: %6.4f\n", BULK_DENS_);
    if ( BULK_DENS_>(0.0334*1.2) ){
      mprintf("\tGIST ATTENTION! water reference density is high, consider using 0.0334 for 1g/cc water density\n");
    }
    if ( BULK_DENS_<(0.0334*0.8) ){
      mprintf("\tGIST: ATTENTION! water reference density is low, consider using 0.0334 for 1g/cc water density\n");
    } 
  }
  else{
//   BULK_DENS_ = 0.033422885325;
    BULK_DENS_ = 0.0334;
    mprintf("\tGIST No water reference density was found, using default: %6.4f, equivalent to 1g/cc\n", BULK_DENS_);
  }

  if ( actionArgs.hasKey("gridcntr") ){
    gridcntr_[0] = actionArgs.getNextDouble(-1);
    gridcntr_[1] = actionArgs.getNextDouble(-1);
    gridcntr_[2] = actionArgs.getNextDouble(-1);
    mprintf("\tGIST grid center: %5.3f %5.3f %5.3f\n", gridcntr_[0],gridcntr_[1],gridcntr_[2]);
  }
  else{
    gridcntr_[0] = 0.0;
    gridcntr_[1] = 0.0;
    gridcntr_[2] = 0.0;
    mprintf("\tGIST No grid center values were found, use default: %5.3f %5.3f %5.3f\n", gridcntr_[0],gridcntr_[1],gridcntr_[2]);
  }

  griddim_.clear();
  griddim_.resize( 3 );
  if ( actionArgs.hasKey("griddim") ){
    griddim_[0] = actionArgs.getNextInteger(-1);
    griddim_[1] = actionArgs.getNextInteger(-1);
    griddim_[2] = actionArgs.getNextInteger(-1);
    mprintf("\tGIST grid dimension: %d %d %d\n", griddim_[0],griddim_[1],griddim_[2]);
  }
  else{
    griddim_[0] = 40;
    griddim_[1] = 40;
    griddim_[2] = 40;
    mprintf("\tGIST No grid dimension values were found, use default: %d %d %d\n", griddim_[0],griddim_[1],griddim_[2]);
  }

  gridspacn_ = actionArgs.getKeyDouble("gridspacn", 0.50);
  mprintf("\tGIST grid spacing: %5.3f A^3\n", gridspacn_);

  mprintf("# Please cite these papers if you use GIST results in a publication: \n");
  mprintf("#   Crystal Nguyen, Michael K. Gilson, and Tom Young,  arXiv:1108.4876v1 (2011)\n");
  mprintf("#   Crystal N. Nguyen, Tom Kurtzman Young, and Michael K. Gilson, J. Chem. Phys. 137, 044101 (2012)\n");
  mprintf("#   Lazaridis, T. Inhomogeneous Fluid Approach to Solvation Thermodynamics. 1. Theory. J. Phys. Chem. B 102, 3531–3541 (1998)\n");

  InitImaging(true);
  
  return Action::OK;
}

// Action_Gist::Setup()
/** Set Gist up for this parmtop. Get masks etc.
  */
Action::RetType Action_Gist::Setup(Topology* currentParm, Topology** parmAddress) {
  CurrentParm_ = currentParm;      
  NFRAME_ = 0;
  max_nwat_ = 0;

  MAX_GRID_PT_ = griddim_[0] * griddim_[1] * griddim_[2];
  Vvox_ = gridspacn_*gridspacn_*gridspacn_;
  G_max_x = griddim_[0] * gridspacn_ + 1.5 ;
  G_max_y = griddim_[1] * gridspacn_ + 1.5 ;
  G_max_z = griddim_[2] * gridspacn_ + 1.5 ;
  
  mprintf("\tGIST number of voxels: %d, voxel volume: %f A^3\n",  MAX_GRID_PT_, Vvox_);

  // Set up grid origin
  gridorig_[0] = gridcntr_[0] - 0.5*griddim_[0]*gridspacn_;
  gridorig_[1] = gridcntr_[1] - 0.5*griddim_[1]*gridspacn_;
  gridorig_[2] = gridcntr_[2] - 0.5*griddim_[2]*gridspacn_;
  mprintf("\tGIST grid origin: %5.3f %5.3f %5.3f\n", gridorig_[0], gridorig_[1], gridorig_[2]);

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
  voxel=0;
  for (int i = 0; i < griddim_[0]; ++i) {
    for (int j = 0; j < griddim_[1]; ++j) {
      for (int k = 0; k < griddim_[2]; ++k) {
        grid_x_[voxel] = Xcrd(i);
        grid_y_[voxel] = Ycrd(j);
        grid_z_[voxel] = Zcrd(k);
        voxel++;
      }
    }
  }

  // wh stands for water-host and can also be written with notation sw
  // ww stands for water-water
  Ewh_dw_.clear();
  Ewh_dw_.resize(MAX_GRID_PT_, 0.0);
  Ewh_norm_.clear();
  Ewh_norm_.resize(MAX_GRID_PT_, 0.0);
  Eww_norm_.clear();
  Eww_norm_.resize(MAX_GRID_PT_, 0.0);
  Eww_dw_.clear();
  Eww_dw_.resize(MAX_GRID_PT_, 0.0);

  if(doEij_) {
    ww_Eij_.clear();
    ww_Eij_.resize(MAX_GRID_PT_);
    for(int i = 1; i < MAX_GRID_PT_; i++) ww_Eij_[i].resize(i);
    
    //CN: need to initialize ww_Eij_ to 0.0 but not Euler angles
    for (int a=1; a<MAX_GRID_PT_; a++)
      for (int l=0; l<a; l++) ww_Eij_[a][l]=0.0;  
  }
  the_vox_.clear();
  the_vox_.resize(MAX_GRID_PT_);
  phi_vox_.clear();
  phi_vox_.resize(MAX_GRID_PT_);
  psi_vox_.clear();
  psi_vox_.resize(MAX_GRID_PT_);

  dTStrans_dw_.clear();
  dTStrans_dw_.resize(MAX_GRID_PT_, 0.0);
  dTStrans_norm_.clear();
  dTStrans_norm_.resize(MAX_GRID_PT_, 0.0); 
  dTSorient_dw_.clear();
  dTSorient_dw_.resize(MAX_GRID_PT_, 0.0);
  dTSorient_norm_.clear();
  dTSorient_norm_.resize(MAX_GRID_PT_, 0.0);

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
  neighbor_dw_.clear();
  neighbor_dw_.resize(MAX_GRID_PT_, 0.0);
  neighbor_norm_.clear();
  neighbor_norm_.resize(MAX_GRID_PT_, 0.0);
  qtet_.clear();
  qtet_.resize(MAX_GRID_PT_, 0.0);
  pol_.clear();
  pol_.resize(MAX_GRID_PT_, 0.0);

  gridwat_.clear();
  gridwat_.resize( currentParm->Nsolvent() );

  // We need box info
  if (currentParm->BoxType() == Box::NOBOX) {
    mprinterr("Error: Gist: Must have explicit solvent with periodic boundaries!");
    return Action::ERR;
  }
  SetupImaging( currentParm->BoxType() );

  resnum =0;
  voxel =0;

  return Action::OK;  
}


// Action_Gist::DoAction()
Action::RetType Action_Gist::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {

  NFRAME_ ++;
//  if (NFRAME_==1) mprintf("GIST Action \n");

  // Simulation box length - assign here because it can vary for npt simulation
  Lx = currentFrame->BoxCrd().BoxX();
  Ly = currentFrame->BoxCrd().BoxY();
  Lz = currentFrame->BoxCrd().BoxZ();
//  if (NFRAME_==1) mprintf("GIST Action box length: %f %f %f \n", Lx, Ly, Lz);
  
  int solventMolecules = CurrentParm_->Nsolvent();
  resnum =0;
  voxel =0;
  resindex1 = 0;
  for (solvmol = CurrentParm_->MolStart();
       solvmol != CurrentParm_->MolEnd(); ++solvmol)
    {
      resindex1++;
      if (!(*solvmol).IsSolvent()) continue;
      Grid( currentFrame );
      voxel = gridwat_[resnum];
      resnum++;   
      NonbondEnergy( currentFrame );
      if (voxel>=MAX_GRID_PT_) continue;
      EulerAngle( currentFrame );
      Dipole( currentFrame );
    }
  if(doOrder_) Order( currentFrame );
  
  //Debugg
//  if (NFRAME_==1) mprintf("GIST  DoAction:  Found %d solvent residues \n", resnum);
  if (solventMolecules != resnum) {
    mprinterr("GIST  DoAction  Error: Number of solvent molecules don't match %d %d\n", solventMolecules, resnum);
  }
  
  return Action::OK;
}

// Action_Gist::NonbondEnergy()
void Action_Gist::NonbondEnergy(Frame *currentFrame) {
  double rij2, rij, r2, r6, r12, f12, f6, e_vdw, e_elec;
  int satom, satom2, atom1, atom2;
  
  int  voxel2 = 0;
  double q1, q2;
  
  // Setup imaging info
  Matrix_3x3 ucell, recip;
  if (ImagingEnabled())
    currentFrame->BoxCrd().ToRecip(ucell, recip);

  // Inner loop has both solute and solvent
  resnum2=0;
  resindex2 = 1;
  // skip if water2 has index larger than water1 so that every pair is only evaluated once
  solvmol2 = CurrentParm_->MolStart();
  for (resindex2=1; resindex2<resindex1; resindex2++)
  {  
    if (!(*solvmol2).IsSolvent()) {
      // Outer loop is not water, break inner loop if water 1 is outside the grid
      if (voxel>=MAX_GRID_PT_) {
        ++solvmol2;
        continue;
      }
    } else { 
      // Inner loop is water
      voxel2 = gridwat_[resnum2];
      resnum2++;
      // skip if both waters are outside the grid
      if (voxel>=MAX_GRID_PT_ && voxel2>=MAX_GRID_PT_) {
        ++solvmol2;
        continue;
      }
    }
    // Loop over all solvent atoms of water 1
    atom1=0;
    for (satom = (*solvmol).BeginAtom(); satom < (*solvmol).EndAtom(); ++satom)
    {
      // Set up coord index for this atom
      const double* XYZ =  currentFrame->XYZ( satom );
  
      atom2=0;
      for (satom2 = (*solvmol2).BeginAtom(); satom2 < (*solvmol2).EndAtom(); ++satom2)
      {    
        // Set up coord index for this atom
        const double* XYZ2 = currentFrame->XYZ( satom2 );
        // Calculate the vector pointing from atom2 to atom1
        rij2 = DIST2(XYZ, XYZ2, ImageType(), currentFrame->BoxCrd(), ucell, recip);
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
        if (!(*solvmol2).IsSolvent()) {
          // solute-solvent interaction
          wh_evdw_[voxel] +=  e_vdw;
          wh_eelec_[voxel] += e_elec;
        } else {
          // solvent-solvent interaction, need to compute for all waters, even 
          // those outside the grid but only one water needs to be inside the grid. 
          if (voxel<MAX_GRID_PT_) {
            ww_evdw_[voxel] +=  e_vdw;
            ww_eelec_[voxel] += e_elec;
            // Store the water neighbor using only O-O distance
            if (atom2==0 && atom1==0 && rij<3.5) {
              neighbor_[voxel] += 1.0;
            }
          }
          // CN: only store Eij[voxel1][voxel2] if both voxels lie on the grid.
          if (voxel2<MAX_GRID_PT_) {
            ww_evdw_[voxel2] +=  e_vdw;
            ww_eelec_[voxel2] += e_elec;
            // Store the water neighbor using only O-O distance
            if (atom2==0 && atom1==0 && rij<3.5) {
              neighbor_[voxel2] += 1.0;
            }
            if(doEij_ && (voxel<MAX_GRID_PT_) ) {
              // full interaction and only print half matrix so if you sums up 
              // all ww_Eij_, you will get the energy of the grid or the system 
              // if the grid covers the entire simulation box. If 1/2 is 
              // included, you will need to multiple the sum by 2. 
              if (voxel>voxel2) {
                ww_Eij_[voxel][voxel2] += e_vdw;
                ww_Eij_[voxel][voxel2] += e_elec;
              } else {
                ww_Eij_[voxel2][voxel] += e_vdw;
                ww_Eij_[voxel2][voxel] += e_elec;
              }
            }//print Eij && voxel<MAX_GRID_PT_
          }
        }//IF is solvent
        atom2++;
      } // END Inner loop ALL atoms
      atom1++;
    } // END Outer loop solvent atoms
    ++solvmol2;
  }  // END Inner loop ALL molecules
}

// Action_Gist::Grid()
void Action_Gist::Grid(Frame *frameIn) {
  int  i, gridindex[3], nH;
  Vec3 comp,  atom_coord;
  i = (*solvmol).BeginAtom();

  gridwat_[resnum] = MAX_GRID_PT_ + 1;
  atom_coord = Vec3(frameIn->XYZ(i));
  // get the components of the water vector
  comp = Vec3(atom_coord) - Vec3(gridorig_);
  nH=0;
  //If Oxygen is far from grid, 1.5A or more in any durection, skip calculation
  if (comp[0]<= G_max_x && comp[1]<= G_max_y && comp[2]<= G_max_z && comp[0]>= -1.5 && comp[1]>= -1.5 && comp[2]>= -1.5 ) {
    //if (comp[0]<= G_max_x || comp[1]<= G_max_y || comp[2]<= G_max_z || comp[0]>= -1.5 || comp[1]>= -1.5 || comp[2]>= -1.5 ) {
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
        voxel = (gridindex[0]*griddim_[1] + gridindex[1])*griddim_[2] + gridindex[2];
        gridwat_[resnum] = voxel;
        nwat_[voxel]++;
        if (max_nwat_ < nwat_[voxel]) max_nwat_ = nwat_[voxel];
      }
    }
    
    // evaluate hydrogen atoms
    for (int a=1; a<=nH; a++) {
      atom_coord = Vec3(frameIn->XYZ(i+a));
      comp = Vec3(atom_coord) - Vec3(gridorig_);
      if (comp[0]<0 || comp[1]<0 || comp[2]<0) continue;
      comp /= gridspacn_;
      gridindex[0] = (int) comp[0];
      gridindex[1] = (int) comp[1];
      gridindex[2] = (int) comp[2];
      if ((gridindex[0]<griddim_[0]) && (gridindex[1]<griddim_[1]) && (gridindex[2]<griddim_[2])) {
        voxel = (gridindex[0]*griddim_[1] + gridindex[1])*griddim_[2] + gridindex[2];
        nH_[voxel]++;
      }
    } 
  }
}

// Action_Gist::EulerAngle()
void Action_Gist::EulerAngle(Frame *frameIn) {

  //if (NFRAME_==1) mprintf("GIST Euler Angles \n");
  Vec3 x_lab, y_lab, z_lab, O_wat, H1_wat, H2_wat, x_wat, y_wat, z_wat, node, v;
  double dp, theta, phi, psi;

  int i = (*solvmol).BeginAtom();
  O_wat = Vec3(frameIn->XYZ(i));
  H1_wat = Vec3(frameIn->XYZ(i+1)) - O_wat;
  H2_wat = Vec3(frameIn->XYZ(i+2)) - O_wat;
  
  // make sure the first three atoms are oxygen followed by two hydrogen
  if ((*CurrentParm_)[i].Element() != Atom::OXYGEN)
    mprintf("Warning: GIST: first coordinates do not belong to oxygen atom\n");
  if ((*CurrentParm_)[i+1].Element() != Atom::HYDROGEN || 
      (*CurrentParm_)[i+2].Element() != Atom::HYDROGEN)
    mprintf("Warning: GIST: second/third coordinates do not belong to hydrogen atoms\n");
  
  // Define lab frame of reference
  x_lab[0]=1.0; x_lab[1]=0; x_lab[2]=0;
  y_lab[0]=0; y_lab[1]=1.0; y_lab[2]=0;
  z_lab[0]=0; z_lab[1]=0; z_lab[2]=1.0;     
  
  // Define the water frame of reference - all axes must be normalized
  // make h1 the water x-axis (but first need to normalized)
  x_wat = H1_wat;
  x_wat.Normalize();
  // the normalized z-axis is the cross product of h1 and h2 
  z_wat = x_wat.Cross( H2_wat );
  z_wat.Normalize();
  // make y-axis as the cross product of h1 and z-axis
  y_wat = z_wat.Cross( x_wat );
  y_wat.Normalize();
  
  // Find the X-convention Z-X'-Z'' Euler angles between the water frame and the lab/host frame
  // First, theta = angle between the water z-axis of the two frames
  dp = z_lab*( z_wat);
  theta = acos(dp);
  //  if (theta>0 && theta<PI) {
  if (theta>1E-5 && theta<Constants::PI-1E-5) {
    // phi = angle between the projection of the water x-axis and the node
    // line of node is where the two xy planes meet = must be perpendicular to both z axes
    // direction of the lines of node = cross product of two normals (z axes)
    // acos of x always gives the angle between 0 and pi, which is okay for theta since theta ranges from 0 to pi
    node = z_lab.Cross( z_wat );
    node.Normalize();
    
    // Second, find the angle phi, which is between x_lab and the node
    dp = node*( x_lab );
    if (dp <= -1.0) phi = Constants::PI;
    else if (dp >= 1.0) phi = Constants::PI;
    else phi = acos(dp);
    // check angle phi
    if (phi>0 && phi<(Constants::TWOPI)) {
      // method 2
      v = x_lab.Cross( node );
      dp = v*( z_lab );
      if (dp<0) phi = Constants::TWOPI - phi;
    }
    
    // Third, rotate the node to x_wat about the z_wat axis by an angle psi
    // psi = angle between x_wat and the node 
    dp = x_wat*( node );
    if (dp<=-1.0) psi = Constants::PI;
    else if (dp>=1.0) psi = 0;
    else psi = acos(dp);
    // check angle psi
    if (psi>0 && psi<(Constants::TWOPI)) {
      // method 2
      Vec3 v = node.Cross( x_wat );
      dp = v*( z_wat );
      if (dp<0) psi = Constants::TWOPI - psi;
    }
    
    if (!(theta<=Constants::PI && theta>=0 && 
          phi<=Constants::TWOPI && phi>=0 && psi<=Constants::TWOPI && psi>=0))
    {
      mprintf("GIST: angles: %f %f %f\n", theta, phi, psi);
      H1_wat.Print("H1_wat");
      H2_wat.Print("H2_wat");
      mprinterr("Error: Euler: angles don't fall into range.\n");
      //break; 
    }
    
    the_vox_[voxel].push_back(theta);
    phi_vox_[voxel].push_back(phi);
    psi_vox_[voxel].push_back(psi);
    nw_angle_[voxel]++;
  }
  //else cout << resnum-1 << " gimbal lock problem, two z_wat paralell" << endl;
} 


// Action_Gist::Dipole()
void Action_Gist::Dipole(Frame *frameIn) {
  
  //if (NFRAME_==1) mprintf("GIST Dipole \n");
  double dipolar_vector[3], charge;
  int satom;
  Vec3 XYZ, sol;

  dipolar_vector[0] = 0.0;
  dipolar_vector[1] = 0.0;
  dipolar_vector[2] = 0.0;
  // Loop over solvent atoms
  for (satom = (*solvmol).BeginAtom(); satom < (*solvmol).EndAtom(); ++satom)
  {
    XYZ =  Vec3(frameIn->XYZ( satom ));
    // Calculate dipole vector. The oxygen of the solvent is used to assign the voxel index to the water.
    // NOTE: the total charge on the solvent should be neutral for this to have any meaning
      
    charge = (*CurrentParm_)[satom].Charge();
    //      cout << resnum-1 << " " << charge << " " << XYZ[0] << " " << XYZ[1] << " " << XYZ[2] << endl;
    dipolar_vector[0] += (charge * XYZ[0]);
    dipolar_vector[1] += (charge * XYZ[1]);
    dipolar_vector[2] += (charge * XYZ[2]);
  }
  dipolex_[voxel] += dipolar_vector[0];
  dipoley_[voxel] += dipolar_vector[1];
  dipolez_[voxel] += dipolar_vector[2];
}

// Action_Gist::Order() 
void Action_Gist::Order(Frame *frameIn) {
//  if (NFRAME_==1) mprintf("GIST Order Parameter \n");
  int i;
  double cos, sum, r1, r2, rij2, x[5], y[5], z[5];
  Vec3 O_wat1, O_wat2, O_wat3, v1, v2;
  resnum=0;

  for (solvmol = CurrentParm_->MolStart();
                              solvmol != CurrentParm_->MolEnd(); ++solvmol)
  {
    if (!(*solvmol).IsSolvent()) continue;

    // obtain 4 closest neighbors for every water
    resnum++;
    voxel = gridwat_[resnum-1];
    if (voxel>=MAX_GRID_PT_) continue;
    // assume that oxygen is the first atom
    i = (*solvmol).BeginAtom();
    O_wat1 = Vec3(frameIn->XYZ( i ));

    r1=1000; r2=1000; resnum2=0;
    for (int a=1; a<5; a++) {
      x[a]=10000;
      y[a]=10000;
      z[a]=10000;
    }
    // Can't make into triangular matrix
    for (solvmol2 = CurrentParm_->MolStart();
                                solvmol2 != CurrentParm_->MolEnd(); ++solvmol2)
    {
      if (!(*solvmol2).IsSolvent()) continue;
      resnum2++;
      if (resnum == resnum2) continue;
      i = (*solvmol2).BeginAtom();
      O_wat2 = Vec3(frameIn->XYZ( i ));      
      rij2 = DIST2_NoImage(O_wat1, O_wat2);
      if (rij2<r1) {
        x[4] = x[3];
        y[4] = y[3];
        z[4] = z[3];
        x[3] = x[2];
        y[3] = y[2];
        z[3] = z[2];
        x[2] = x[1];
        y[2] = y[1];
        z[2] = z[1];
        r1 = rij2;
        x[1] = O_wat2[0];
        y[1] = O_wat2[1];
        z[1] = O_wat2[2];
      }
    }
    
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
    qtet_[voxel] += (1.0 - (3.0/8)*sum);
  }
}


void Action_Gist::Print() {
  
  // Implement NN to compute orientational entropy for each voxel
  double NNr, rx, ry, rz, rR, dbl;
  dTSorienttot_=0;
  for (int gr_pt=0; gr_pt<MAX_GRID_PT_; gr_pt++) {
    dTSorient_norm_[gr_pt]=0;  
    int nwtot = nw_angle_[gr_pt];
    if (nwtot<=1) continue;
    for (int n=0; n<nwtot; n++) {
      NNr=10000;
      for (int l=0; l<nwtot; l++) {
        if (l==n) continue;
        rx = cos(the_vox_[gr_pt][l]) - cos(the_vox_[gr_pt][n]);
        ry = phi_vox_[gr_pt][l] - phi_vox_[gr_pt][n];
        rz = psi_vox_[gr_pt][l] - psi_vox_[gr_pt][n];
        if (ry>Constants::PI) ry = Constants::TWOPI-ry;
        else if (ry<-Constants::PI) ry = Constants::TWOPI+ry;
        if (rz>Constants::PI) rz = Constants::TWOPI-rz;
        else if (rz<-Constants::PI) rz = Constants::TWOPI+rz;
        rR = sqrt(rx*rx + ry*ry + rz*rz);
        if (rR>0 && rR<NNr) NNr = rR;
      }
      if (NNr<9999 && NNr>0) {
        dbl = log(NNr*NNr*NNr*nwtot/(3.0*Constants::TWOPI));
        dTSorient_norm_[gr_pt] += dbl;
      }
    }
    dTSorient_norm_[gr_pt] = Constants::GASK_KCAL*300*(dTSorient_norm_[gr_pt]/nwtot+0.5772);
    dTSorient_dw_[gr_pt] = dTSorient_norm_[gr_pt]*nwat_[gr_pt]/(NFRAME_*Vvox_);
    dTSorienttot_ += dTSorient_dw_[gr_pt];
  }
  dTSorienttot_ *= Vvox_;
  mprintf("Maximum number of waters found in one voxel for %d frames = %d\n", NFRAME_, max_nwat_);
  mprintf("Total referenced orientational entropy of the grid: dTS_orient = %9.5f kcal/mol\n", dTSorienttot_);
  
  // Compute translational entropy for each voxel
  dTStranstot_=0;
  for (int a=0; a<MAX_GRID_PT_; a++) {
    dens_[a] = 1.0*nwat_[a]/(NFRAME_*Vvox_);
    g_[a] = dens_[a]/BULK_DENS_;
    gH_[a] = 1.0*nH_[a]/(NFRAME_*Vvox_*2*BULK_DENS_);
    if (nwat_[a]>1) {
       dTStrans_dw_[a] = -Constants::GASK_KCAL*BULK_DENS_*300*g_[a]*log(g_[a]);
       // unreferenced translational entropy
       //TStrans_dw_[a] = -Constants::GASK_KCAL*300*dens_[a]*log(dens_[a]);
       dTStrans_norm_[a] = dTStrans_dw_[a]/dens_[a];
       dTStranstot_ += dTStrans_dw_[a];
    }
    else {
       dTStrans_dw_[a]=0; dTStrans_norm_[a]=0;
    }
  }
  dTStranstot_ *= Vvox_;
  mprintf("Total referenced translational entropy of the grid: dTS_trans = %9.5f kcal/mol\n", dTStranstot_);

  // Compute average unreferenced energy
  for (int a=0; a<MAX_GRID_PT_; a++) {
    if (nwat_[a]>0) {
       Ewh_dw_[a] = (wh_evdw_[a]+wh_eelec_[a])/(NFRAME_*Vvox_);
       Ewh_norm_[a] = (wh_evdw_[a]+wh_eelec_[a])/nwat_[a];
       Eww_dw_[a] = (ww_evdw_[a]+ww_eelec_[a])/(2*NFRAME_*Vvox_);
       Eww_norm_[a] = (ww_evdw_[a]+ww_eelec_[a])/(2*nwat_[a]); 
       Ewhtot_ += Ewh_dw_[a];
       Ewwtot_ += Eww_dw_[a];
    }
    else {
       Ewh_dw_[a]=0; Ewh_norm_[a]=0; Eww_norm_[a]=0; Eww_dw_[a]=0;
    }

    // Compute the average number of water neighbor, average order parameter, and average dipole density 
    if (nwat_[a]>0) {
      qtet_[a] /= nwat_[a];
      neighbor_norm_[a] = 1.0*neighbor_[a]/nwat_[a];
    }
    neighbor_dw_[a] = 1.0*neighbor_[a]/(NFRAME_*Vvox_);

    dipolex_[a] /= (0.20822678*NFRAME_*Vvox_);
    dipoley_[a] /= (0.20822678*NFRAME_*Vvox_);
    dipolez_[a] /= (0.20822678*NFRAME_*Vvox_);

    pol_[a] = sqrt(dipolex_[a]*dipolex_[a] + dipoley_[a]*dipoley_[a] + dipolez_[a]*dipolez_[a]);    
  }
  Ewhtot_ *= Vvox_;
  Ewwtot_ *= Vvox_;
  mprintf("Total water-host energy of the grid: Esw = %9.5f kcal/mol\n", Ewhtot_);
  mprintf("Total unreferenced water-water energy of the grid: Eww = %9.5f kcal/mol\n", Ewwtot_);  

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
  PrintDX("gist-Esw-dens.dx", Ewh_dw_);
  PrintDX("gist-Eww-dens.dx", Eww_dw_);
  PrintDX("gist-dTStrans-dens.dx", dTStrans_dw_);
  PrintDX("gist-dTSorient-dens.dx", dTSorient_dw_);
  PrintDX("gist-neighbor-norm.dx", neighbor_norm_); 
  PrintDX("gist-dipole-dens.dx", pol_);
  //stored as doubles
  PrintDX_double("gist-order-norm.dx", qtet_);
  PrintDX_double("gist-dipolex-dens.dx", dipolex_);
  PrintDX_double("gist-dipoley-dens.dx", dipoley_);
  PrintDX_double("gist-dipolez-dens.dx", dipolez_);
  

  if (!datafile_.empty()) {
    PrintOutput(datafile_);
  }
  else{
    PrintOutput("gist-output.dat");
  }
//  gist_t_end=clock();
//  gist_t_diff = diffclock(gist_t_end,gist_t_begin);
//  cout << "GIST Time elapsed: " << gist_t_diff << " ms"<< endl;
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
  outfile.Printf("GIST Output, information printed per voxel\n");
  outfile.Printf("voxel xcoord ycoord zcoord population g_O g_H ");
  outfile.Printf("dTStrans-dens(kcal/mol/A^3) dTStrans-norm(kcal/mol) dTSorient-dens(kcal/mol/A^3) dTSorient-norm(kcal/mol) ");
  outfile.Printf("Esw-dens(kcal/mol/A^3) Esw-norm(kcal/mol) ");
  outfile.Printf("Eww-dens(kcal/mol/A^3) Eww-norm(kcal/mol) ");
  outfile.Printf("Dipole_x-dens(D/A^3) Dipole_y-dens(D/A^3) Dipole_z-dens(D/A^3) Dipole-dens(D/A^3) neighbor-dens(1/A^3) neighbor-norm order-norm\n");
  // Now print out the data. 
  for (int i=0; i<MAX_GRID_PT_; i++){
    outfile.Printf( "%d %g %g %g %d %g %g ",i , grid_x_[i] , grid_y_[i], grid_z_[i], nwat_[i] , g_[i], gH_[i] );
    outfile.Printf( "%g %g %g %g ",dTStrans_dw_[i], dTStrans_norm_[i], dTSorient_dw_[i] , dTSorient_norm_[i]);
    outfile.Printf( "%g %g ",Ewh_dw_[i], Ewh_norm_[i] );
    outfile.Printf( "%g %g ",Eww_dw_[i] , Eww_norm_[i] );
    outfile.Printf( "%g %g %g %g ",dipolex_[i] , dipoley_[i] , dipolez_[i] , pol_[i] );
    outfile.Printf( "%g %g %g \n",neighbor_dw_[i] , neighbor_norm_[i] , qtet_[i]);
    }
  outfile.CloseFile();

  if(doEij_){
    // print half interaction to be consistent with Eww term. However, only print half matrix (triangular instead of square). Thus, to get the interaction of the grid, must multiply the sum by 2. 
    if (outfile.OpenWrite("Eww_ij.dat") == 0) {
      double dbl, Eijtot=0;
        for (int a=1; a<MAX_GRID_PT_; a++) {
        for (int l=0; l<a; l++) {
          dbl = ww_Eij_[a][l];
          if (dbl!=0) {
            dbl /= (2*NFRAME_);
            outfile.Printf("%10d %10d %12.5E\n", a, l, dbl);
            Eijtot += dbl;
          }
        }
      }
      outfile.Printf("Eijtot %12.5E\n", Eijtot);
      outfile.CloseFile();
    } else
      mprinterr("Error: Could not open 'Eww_Eij.dat' for writing.\n");
  }
}
