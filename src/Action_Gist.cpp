// Gist 
#include <cmath>
#include <iostream> // cout
#include <cstring> // memset
#include <fstream>
using namespace std;
#include "Action_Gist.h"
#include "CpptrajFile.h"
#include "CpptrajStdio.h"
#include "Constants.h" // RADDEG
#include "DistRoutines.h"
#include "DataSet_integer.h"
#include "Box.h"
#include "StringRoutines.h" 
#include "Constants.h" // GASCONSTANT

// CONSTRUCTOR
Action_Gist::Action_Gist() :
  CurrentParm_(0),
  kes_(1.0),
  ELJ_(0),
  Eelec_(0),
  watermodel_(false),
  useTIP3P_(false),
  useTIP4P_(false),
  useTIP4PEW_(false),
  griddim_(0)
//  gridindex_(0)
{
  mprintf("\tGIST: INIT \n");
  gridcntr_[0] = -1;
  gridcntr_[1] = -1;
  gridcntr_[2] = -1;
    
  gridorig_[0] = -1;
  gridorig_[1] = -1;
  gridorig_[2] = -1;
  
  gridspacn_ = 0;
 } 


void Action_Gist::Help() {
  mprintf("gist <watermodel>[{tip3p|tip4p|tip4pew}] [gridcntr <xval> <yval> <zval>] [griddim <xval> <yval> <zval>] [gridspacn <spaceval>] [out <filename>] \n");
  mprintf("\tCalculate GIST between water molecules in selected site \n");
}


// Action_Gist::init()
Action::RetType Action_Gist::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
				  DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  mprintf("\tGIST: init \n");
  // Get keywords
  
  // Dataset to store gist results
  datafile_ = actionArgs.GetStringKey("out");
  // Generate the data set name, and hold onto the master data set list
  std::string ds_name = actionArgs.GetStringKey("name");
  ds_name = myDSL_.GenerateDefaultName("GIST");
  // We have 4?? data sets Add them here
  // Now add all of the data sets
  for (int i = 0; i < 4; i++) {
    myDSL_.AddSetAspect(DataSet::DOUBLE, ds_name,
			integerToString(i+1).c_str());
  }
  //  myDSL_.AddSet(DataSet::DOUBLE, ds_name, NULL);
  
  useTIP3P_ = actionArgs.hasKey("tip3p");
  useTIP4P_ = actionArgs.hasKey("tip4p");
  useTIP4PEW_ = actionArgs.hasKey("tip4pew");
  if (!useTIP3P_ && !useTIP4P_ && !useTIP4PEW_) {
    mprinterr("Error: gist: Only water models supprted are TIP3P and TIP4P\n");
    return Action::ERR;
  }

  // Set Bulk Energy based on water model
  if (useTIP3P_) BULK_E_ = -19.0653;
  if (useTIP4P_ || useTIP4PEW_) BULK_E_ = -22.06;
  
  // Set Bulk Density
  BULK_DENS_ = 0.033422885325;

  if ( actionArgs.hasKey("gridcntr") ){
    gridcntr_[0] = actionArgs.getNextDouble(-1);
    gridcntr_[1] = actionArgs.getNextDouble(-1);
    gridcntr_[2] = actionArgs.getNextDouble(-1);
    mprintf("\tGIST grid center: %5.3f %5.3f %5.3f\n", gridcntr_[0],gridcntr_[1],gridcntr_[2]);
  }
  else{
    mprintf("\tGIST: No grid center values were found, using default\n");
    gridcntr_[0] = 0.0;
    gridcntr_[1] = 0.0;
    gridcntr_[2] = 0.0;
    mprintf("\tGIST grid center: %5.3f %5.3f %5.3f\n", gridcntr_[0],gridcntr_[1],gridcntr_[2]);
  }

  griddim_ = new int[3];
  memset(griddim_, 0, 3 * sizeof(int));
  if ( actionArgs.hasKey("griddim") ){
    griddim_[0] = actionArgs.getNextInteger(-1);
    griddim_[1] = actionArgs.getNextInteger(-1);
    griddim_[2] = actionArgs.getNextInteger(-1);
    mprintf("\tGIST grid dimension: %d %d %d \n", griddim_[0],griddim_[1],griddim_[2]);
  }
  else{
    mprintf("\tGIST: No grid dimensiom values were found, using default (box size) \n");
    griddim_[0] = 100;
    griddim_[1] = 100;
    griddim_[2] = 100;
    mprintf("\tGIST grid dimension: %d %d %d \n", griddim_[0],griddim_[1],griddim_[2]);
  }


  gridspacn_ = actionArgs.getKeyDouble("gridspacn", 0.50);
  mprintf("\tGIST grid spacing: %5.3f \n", gridspacn_);

  return Action::OK;
}

// Action_Gist::setup()
/** Set Gist up for this parmtop. Get masks etc.
  */
Action::RetType Action_Gist::Setup(Topology* currentParm, Topology** parmAddress) {
  mprintf("GIST Setup \n");

  CurrentParm_ = currentParm;      
  NFRAME_ = 0;
  max_nwat_ = 0;

  MAX_GRID_PT_ = griddim_[0] * griddim_[1] * griddim_[2];
  Vvox_ = gridspacn_*gridspacn_*gridspacn_;
  mprintf("\tGIST Setup: %d %d %d %d %f \n", griddim_[0], griddim_[1], griddim_[2], MAX_GRID_PT_, Vvox_);

  // Set up cumulative energy arrays
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

  dEwh_dw_.clear();
  dEwh_dw_.resize(MAX_GRID_PT_, 0.0);
  dEww_dw_ref_.clear();
  dEww_dw_ref_.resize(MAX_GRID_PT_, 0.0);
  dEwh_norm_.clear();
  dEwh_norm_.resize(MAX_GRID_PT_, 0.0);
  dEww_norm_ref_.clear();
  dEww_norm_ref_.resize(MAX_GRID_PT_, 0.0);

  ww_Eij_.clear();
  ww_Eij_.resize(MAX_GRID_PT_);
  for(int i = 0; i < MAX_GRID_PT_; i++) ww_Eij_[i].resize(MAX_GRID_PT_);

  //CN: need to initialize ww_Eij_ to 0.0 but not Euler angles
  for (int a=0; a<MAX_GRID_PT_; a++)
    for (int l=0; l<MAX_GRID_PT_; l++) ww_Eij_[a][l]=0.0;
  
  the_vox_.clear();
  the_vox_.resize(MAX_GRID_PT_);
  phi_vox_.clear();
  phi_vox_.resize(MAX_GRID_PT_);
  psi_vox_.clear();
  psi_vox_.resize(MAX_GRID_PT_);

  TStrans_dw_.clear();
  TStrans_dw_.resize(MAX_GRID_PT_, 0.0);
  TStrans_norm_.clear();
  TStrans_norm_.resize(MAX_GRID_PT_, 0.0); 
  TSNN_.clear();
  TSNN_.resize(MAX_GRID_PT_, 0.0);
  TSwNN_.clear();
  TSwNN_.resize(MAX_GRID_PT_, 0.0);

  nwat_.clear();
  nwat_.resize(MAX_GRID_PT_, 0);
  nw_angle_.clear();
  nw_angle_.resize(MAX_GRID_PT_, 0);
  dens_.clear();
  dens_.resize(MAX_GRID_PT_, 0);
  g_.clear();
  g_.resize(MAX_GRID_PT_, 0);

  //Set atom charges
  atom_charge_.clear();
  atom_charge_.reserve( currentParm->Natom() );
  for (Topology::atom_iterator atom = currentParm->begin(); atom != currentParm->end(); ++atom)
    atom_charge_.push_back( (*atom).Charge() * ELECTOAMBER );
  gridwat_.clear();
  gridwat_.reserve( currentParm->Natom() );


  // Set up grid origin
  gridorig_[0] = gridcntr_[0] - 0.5*griddim_[0]*gridspacn_;
  gridorig_[1] = gridcntr_[1] - 0.5*griddim_[1]*gridspacn_;
  gridorig_[2] = gridcntr_[2] - 0.5*griddim_[2]*gridspacn_;
  mprintf("\tGIST grid origin: %5.3f %5.3f %5.3f\n", gridorig_[0],gridorig_[1],gridorig_[2]);
 
  return Action::OK;  
}


// Action_Gist::action()
Action::RetType Action_Gist::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {

  if (NFRAME_==1) mprintf("GIST Action \n");

  //select water molecules
  NFRAME_++;
  Grid( currentFrame, CurrentParm_);
  NonbondEnergy( currentFrame, CurrentParm_ );
  EulerAngle( currentFrame, CurrentParm_);

  return Action::OK;
}


static void GetLJparam(Topology const& top, double& A, double& B, 
                              int atom1, int atom2)
{
  // In Cpptraj, atom numbers start from 1, so subtract 1 from the NB index array
  int param = (top.Ntypes() * (top[atom1].TypeIndex()-1)) + top[atom2].TypeIndex()-1;
  int index = top.NB_index()[param] - 1;
  A = top.LJA()[index];
  B = top.LJB()[index];
}

void Action_Gist::NonbondEnergy(Frame *currentFrame, Topology *CurrentParm) {
  double delta2, Acoef, Bcoef, deltatest;
  Vec3 XYZ,XYZ2, JI;
  double rij2,rij,r2,r6,r12,f12,f6,e_vdw,qiqj,e_elec;
  int satom,satom2;
  
  Topology::mol_iterator solvmol, solvmol2;
  if (NFRAME_==1) mprintf("GIST NonbondEnergy \n");
  int voxel,voxel2;
  resnum=0;
  ELJ_ = 0.0;
  Eelec_ = 0.0;
  kes_ = 0.0;

  // Loop over all molecules
  // Outer loop
  for (solvmol = CurrentParm_->MolStart();
       solvmol != CurrentParm_->MolEnd(); ++solvmol)
    {
      if (!(*solvmol).IsSolvent()) continue;
      voxel = gridwat_[resnum];
      // if main water is outside the grid, skip
      resnum++;
      if (voxel>=MAX_GRID_PT_) continue;
      // Loop over solvent atoms
      for (satom = (*solvmol).BeginAtom(); satom < (*solvmol).EndAtom(); ++satom)
	{
	  // Set up coord index for this atom
	  XYZ =  Vec3(currentFrame->XYZ( satom ));	  
	  // Inner loop	  
	  resnum2=0;
	  for (solvmol2 = CurrentParm_->MolStart();
	       solvmol2 != CurrentParm_->MolEnd(); ++solvmol2)
	    {
	      if (!(*solvmol2).IsSolvent()) { //notsolvent loop
		for (satom2 = (*solvmol2).BeginAtom(); satom2 < (*solvmol2).EndAtom(); ++satom2)
		  {
		    
		    // Set up coord index for this atom
		    XYZ2 = Vec3(currentFrame->XYZ( satom2 ));
		    // Calculate the vector pointing from atom2 to atom1
		    JI = Vec3(XYZ) - Vec3(XYZ2);
		    rij2 = JI.Magnitude2();
		    // Normalize
		    rij = sqrt(rij2);
		    JI /= rij;
		    // LJ energy 
		    GetLJparam(*CurrentParm, Acoef, Bcoef, satom, satom2);
		    r2    = 1 / rij2;
		    r6    = r2 * r2 * r2;
		    r12   = r6 * r6;
		    f12   = Acoef * r12;  // A/r^12
		    f6    = Bcoef * r6;   // B/r^6
		    e_vdw = f12 - f6;     // (A/r^12)-(B/r^6)
		    ELJ_ = ELJ_ + e_vdw;
		    // LJ Force 
		    //force=((12*f12)-(6*f6))*r2; // (12A/r^13)-(6B/r^7)
		    // Coulomb energy 
		    qiqj = atom_charge_[satom] * atom_charge_[satom2];
		    e_elec = kes_ * (qiqj/rij);
		    Eelec_ = Eelec_ + e_elec;
		    // Cumulative evdw - divide between both atoms
//		    delta2 = e_vdw * 0.5;
		    wh_evdw_[voxel] += e_vdw;
//		    deltatest = delta2;
		    // Cumulative eelec - divide between both atoms
//		    delta2 = e_elec * 0.5;
		    wh_eelec_[voxel] += e_elec;
		    //mprintf("GIST Action notsolvent atom1 %d atom2 %d res1 %d res2 %d vdW %f eelec %f \n",satom,satom2,resnum-1,resnum2,e_vdw,e_elec);
		    
		    // ----------------------------------------
		  } // END Inner loop non-solvent atoms
	      } //If loop notsolvent loop
	      else{ // Solvent loop
		// skip water intearction with with itself
		voxel2 = gridwat_[resnum2];
		resnum2++;
		if ((resnum-1)==(resnum2-1)) continue;
		for (satom2 = (*solvmol2).BeginAtom(); satom2 < (*solvmol2).EndAtom(); ++satom2)
		  {
		    
		    // Set up coord index for this atom
		    XYZ2 = Vec3(currentFrame->XYZ( satom2 ));
		    // Calculate the vector pointing from atom2 to atom1
		    JI = Vec3(XYZ) - Vec3(XYZ2);
		    rij2 = JI.Magnitude2();
		    // Normalize
		    rij = sqrt(rij2);
		    JI /= rij;
		    // LJ energy 
		    GetLJparam(*CurrentParm, Acoef, Bcoef, satom, satom2);
		    r2    = 1 / rij2;
		    r6    = r2 * r2 * r2;
		    r12   = r6 * r6;
		    f12   = Acoef * r12;  // A/r^12
		    f6    = Bcoef * r6;   // B/r^6
		    e_vdw = f12 - f6;     // (A/r^12)-(B/r^6)
		    //ELJ_ = ELJ_ + e_vdw;
		    // LJ Force 
		    //force=((12*f12)-(6*f6))*r2; // (12A/r^13)-(6B/r^7)
		    //scalarmult(f,JI,F);
		    // Coulomb energy 
		    qiqj = atom_charge_[satom] * atom_charge_[satom2];
		    e_elec = kes_ * (qiqj/rij);
		    //Eelec_ = Eelec_ + e_elec;
		    // Cumulative evdw - divide between both atoms
		    // CN: Here we divide each energy pair by 1/2 because each pairs is counted twice. So it's still the full interactino. 
		    //delta2 = e_vdw * 0.5;
		    //ww_evdw_[voxel] +=  delta2;
		    ww_evdw_[voxel] +=  e_vdw;
		    // CN: only store Eij[voxel1][voxel2] if both voxels lie on the grid.
		    if (voxel2<MAX_GRID_PT_) ww_Eij_[voxel][voxel2] += delta2;
		    // Cumulative eelec - divide between both atoms
		    //delta2 = e_elec * 0.5;
		    //ww_eelec_[voxel] += delta2;
		    ww_eelec_[voxel] += e_elec;
		    if (voxel2<MAX_GRID_PT_) ww_Eij_[voxel][voxel2] += delta2;
		    //mprintf("GIST Action solvent atom1 %d atom2 %d res1 %d res2 %d vdW %f eelec %f \n",satom,satom2,resnum-1,resnum2-1,e_vdw,e_elec);
		    
		    // ----------------------------------------
		  } // END Inner loop solvent atoms
	      } //Else solvent loop
	    } // END Inner loop ALL molecules 
	} // END Outer loop solvent atoms
    } // END Outer loop solvent molecules
}


// Action_Gist::Grid()
void Action_Gist::Grid(Frame *frameIn, Topology* CurrentParm) {

  resnum=0;
  int voxel, i, gridindex[3];
  Vec3 comp, compnew, O_wat;
  double rij;

  for (Topology::mol_iterator solvmol = CurrentParm_->MolStart();
       solvmol != CurrentParm_->MolEnd(); ++solvmol)
    {
      if (!(*solvmol).IsSolvent()) continue;
      i = (*solvmol).BeginAtom();
      O_wat = Vec3(frameIn->XYZ(i));
      // get the components of the water vector
      comp = Vec3(O_wat) - Vec3(gridorig_);
      rij = sqrt( comp.Magnitude2() );
      compnew = comp/gridspacn_;
      gridindex[0] = (int) compnew[0];
      gridindex[1] = (int) compnew[1];
      gridindex[2] = (int) compnew[2];
      gridwat_[resnum] = 100000000;
      if (gridindex[0]>=0 && gridindex[1]>=0 && gridindex[2]>=0 && (gridindex[0]<griddim_[0]) && (gridindex[1]<griddim_[1]) && (gridindex[2]<griddim_[2]))
        {
	  // this water belongs to grid point gridindex[0], gridindex[1], gridindex[2]
          voxel = (gridindex[0]*griddim_[1] + gridindex[1])*griddim_[2] + gridindex[2];
	  //voxel = (((gridindex_[2]-1)*griddim_[1] + gridindex_[1])-1)*griddim_[0] + gridindex_[0];
          gridwat_[resnum] = voxel;
          //mprintf("fm=%d, resnum=%d, voxel=%d %d\n", NFRAME_, resnum, voxel, gridwat_[resnum]);
        }
      resnum++;
    }

    //Debugg
    int solventMolecules = CurrentParm_->Nsolvent();
    if (NFRAME_==1) mprintf("GIST  Grid:  Found %d solvent residues \n", resnum);
    if (solventMolecules != resnum) {
      mprinterr("GIST  Grid  Error: No solvent molecules don't match %d %d\n", solventMolecules, resnum);
    }
}
	

void Action_Gist::EulerAngle(Frame *frameIn, Topology* CurrentParm) {

  Vec3 O_wat, H1_wat, H2_wat, x_wat, y_wat, z_wat, node, v;
  int voxel, resnum3=0;
  double cp, dp, theta, phi, psi;
  for (Topology::mol_iterator solvmol = CurrentParm_->MolStart();
       solvmol != CurrentParm_->MolEnd(); ++solvmol) {
    // if molecule is not a solvent, skip
    if (!(*solvmol).IsSolvent()) continue;
    voxel = gridwat_[resnum3];
    resnum3++;
    if (voxel >= MAX_GRID_PT_) continue;
    nwat_[voxel] ++;
    if (max_nwat_ < nwat_[voxel]) {
      max_nwat_ = nwat_[voxel];
    }
      
    int i = (*solvmol).BeginAtom();
    O_wat = Vec3(frameIn->XYZ(i));
    H1_wat = Vec3(frameIn->XYZ(i+1)) - O_wat;
    H2_wat = Vec3(frameIn->XYZ(i+2)) - O_wat;

    // Define lab frame of reference
    x_lab[0]=1.0; x_lab[1]=0; x_lab[2]=0;
    y_lab[0]=0; y_lab[1]=1.0; y_lab[2]=0;
    z_lab[0]=0; z_lab[1]=0; z_lab[2]=1.0;     
      
    // Define the water frame of reference - all axes must be normalized
    // make h1 the water x-axis (but first need to normalized)
    x_wat = H1_wat;
    cp = x_wat.Normalize();
    // the normalized z-axis is the cross product of h1 and h2 
    z_wat = x_wat.Cross( H2_wat );
    cp = z_wat.Normalize();
    // make y-axis as the cross product of h1 and z-axis
    y_wat = z_wat.Cross( x_wat );
    cp = y_wat.Normalize();
      
    // Find the X-convention Z-X'-Z'' Euler angles between the water frame and the lab/host frame
    // First, theta = angle between the water z-axis of the two frames
    dp = z_lab*( z_wat);
    theta = acos(dp);
    if (theta>0 && theta<PI) {
      // phi = angle between the projection of the water x-axis and the node
      // line of node is where the two xy planes meet = must be perpendicular to both z axes
      // direction of the lines of node = cross product of two normals (z axes)
      // acos of x always gives the angle between 0 and pi, which is okay for theta since theta ranges from 0 to pi
      node = z_lab.Cross( z_wat );
      cp = node.Normalize();

      // Second, find the angle phi, which is between x_lab and the node
      dp = node*( x_lab );
      if (dp <= -1.0) phi = PI;
      else if (dp >= 1.0) phi = PI;
      else phi = acos(dp);
      // check angle phi
      if (phi>0 && phi<(2*PI)) {
        // method 2
	v = x_lab.Cross( node );
	dp = v*( z_lab );
	if (dp<0) phi = 2*PI - phi;
      }

      // Third, rotate the node to x_wat about the z_wat axis by an angle psi
      // psi = angle between x_wat and the node 
      dp = x_wat*( node );
      if (dp<=-1.0) psi = PI;
      else if (dp>=1.0) psi = 0;
      else psi = acos(dp);
      // check angle psi
      if (psi>0 && psi<(2*PI)) {
        // method 2
        Vec3 v = node.Cross( x_wat );
        dp = v*( z_wat );
	if (dp<0) psi = 2*PI - psi;
      }
	
      // DEBUG
      // The total rotational matrix for transforming the water frame onto the lab frame
      float ** mat_W = new float * [3];
      for (int a=0; a<3; a++) {
	mat_W[a] = new float [3];
      }
      mat_W[0][0] = cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi);
      mat_W[0][1] = cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi);
      mat_W[0][2] = sin(psi)*sin(theta);
      mat_W[1][0] = -sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi);
      mat_W[1][1] = -sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi);
      mat_W[1][2] = cos(psi)*sin(theta);
      mat_W[2][0] = sin(theta)*sin(phi);
      mat_W[2][1] = -sin(theta)*cos(phi);
      mat_W[2][2] = cos(theta);

      // apply the rotational matrix to the water frame of reference
//      Vec3 x_res = x_wat*( mat_W );
//      Vec3 y_res = y_wat*( mat_W );
//      Vec3 z_res = z_wat*( mat_W );
      //I uncommented this since it won't take a matrix * Vec3 operation and this is equivalent, right?
      x_res[0] = x_wat[0]*mat_W[0][0] + x_wat[1]*mat_W[0][1] + x_wat[2]*mat_W[0][2];
      x_res[1] = x_wat[0]*mat_W[1][0] + x_wat[1]*mat_W[1][1] + x_wat[2]*mat_W[1][2];
      x_res[2] = x_wat[0]*mat_W[2][0] + x_wat[1]*mat_W[2][1] + x_wat[2]*mat_W[2][2];
      y_res[0] = y_wat[0]*mat_W[0][0] + y_wat[1]*mat_W[0][1] + y_wat[2]*mat_W[0][2];
      y_res[1] = y_wat[0]*mat_W[1][0] + y_wat[1]*mat_W[1][1] + y_wat[2]*mat_W[1][2];
      y_res[2] = y_wat[0]*mat_W[2][0] + y_wat[1]*mat_W[2][1] + y_wat[2]*mat_W[2][2];
      z_res[0] = z_wat[0]*mat_W[0][0] + z_wat[1]*mat_W[0][1] + z_wat[2]*mat_W[0][2];
      z_res[1] = z_wat[0]*mat_W[1][0] + z_wat[1]*mat_W[1][1] + z_wat[2]*mat_W[1][2];
      z_res[2] = z_wat[0]*mat_W[2][0] + z_wat[1]*mat_W[2][1] + z_wat[2]*mat_W[2][2];

      for (int a=0; a<3; a++) {
	delete [] mat_W[a];
      }
      delete [] mat_W;

      double rRx = x_res*( x_lab );
      double rRy = y_res*( y_lab );
      double rRz = z_res*( z_lab );
      if (rRx>1+1E-6 || rRx<1-1E-6 || rRy>1+1E-6 || rRy<1-1E-6 || rRz>1+1E-6 || rRz<1-1E-6) {
        std::cout  << "wat=" << resnum3-1 << ", gr=" << voxel << " ROTATION IS BAD!" << std::endl;
        std::cout << "rx=" << rRx << ", ry=" << rRy << ", rz=" << rRz << std::endl;
        std::cout << "water new x axis: " << x_res[0] << " " << x_res[1] << " " << x_res[2] << std::endl;
        std::cout << "water new y axis: " << y_res[0] << " " << y_res[1] << " " << y_res[2] << std::endl;
        std::cout << "water new z axis: " << z_res[0] << " " << z_res[1] << " " << z_res[2] << std::endl;
        mprinterr("Error: Euler: BAD ROTATION.\n");
        break;	
      }
	
      if (!(theta<=PI && theta>=0 && phi<=2*PI && phi>=0 && psi<=2*PI && psi>=0)) {
        std::cout << "angles: " << theta << " " << phi << " " << psi << std::endl;
        std::cout << H1_wat[0] << " " << H1_wat[1] << " " << H1_wat[2] << " " << H2_wat[0] << " " << H2_wat[1] << " " << H2_wat[2] << std::endl;
        mprinterr("Error: Euler: angles don't fall into range.\n");
        break; 
      }
   
      the_vox_[voxel].push_back(theta);
      phi_vox_[voxel].push_back(phi);
      psi_vox_[voxel].push_back(psi);
      nw_angle_[voxel]++;
    }
    else std::cout << " " << resnum3-1 << " gimbal lock problem, two z_wat paralell" << std::endl;
  }
} 


void Action_Gist::Print() {

  // Implement NN to compute orientational entropy for each voxel
  double NNr, rx, ry, rz, rR, dbl;
  TSNNtot_=0;
  for (int gr_pt=0; gr_pt<MAX_GRID_PT_; gr_pt++) {
    TSNN_[gr_pt]=0; TSwNN_[gr_pt]=0;  
    int nwtot = nw_angle_[gr_pt];
    if (nwtot<=1) continue;
    for (int n=0; n<nwtot; n++) {
      NNr=10000;
      for (int l=0; l<nwtot; l++) {
        if (l==n) continue;
        rx = cos(the_vox_[gr_pt][l]) - cos(the_vox_[gr_pt][n]);
        ry = phi_vox_[gr_pt][l] - phi_vox_[gr_pt][n];
        rz = psi_vox_[gr_pt][l] - psi_vox_[gr_pt][n];
        if (ry>PI) ry = 2*PI-ry;
        else if (ry<-PI) ry = 2*PI+ry;
        if (rz>PI) rz = 2*PI-rz;
        else if (rz<-PI) rz = 2*PI+rz;
        rR = sqrt(rx*rx + ry*ry + rz*rz);
        if (rR>0 && rR<NNr) NNr = rR;
      }
      if (NNr<9999 && NNr>0) {
	dbl = log(NNr*NNr*NNr*nwtot/(3*2*PI));
	TSwNN_[gr_pt] += dbl;
      }	
    }
    TSwNN_[gr_pt] = GASCONSTANT*0.3*0.239*(TSwNN_[gr_pt]/nwtot+0.5772);
    TSNN_[gr_pt] = TSwNN_[gr_pt]*nwat_[gr_pt]/(NFRAME_*Vvox_);
    TSNNtot_ += TSNN_[gr_pt];
  }
  TSNNtot_ *= Vvox_;
  mprintf("Max number of water = %d \n", max_nwat_);
  mprintf("Total orientational entropy of the grid: S_NN = %9.5f kcal/mol, Nf=%d\n", TSNNtot_, NFRAME_);
  
  // Compute translational entropy for each voxel
  TStranstot_=0;
  for (int a=0; a<MAX_GRID_PT_; a++) {
    dens_[a] = nwat_[a]/(NFRAME_*Vvox_);
    g_[a] = dens_[a]/BULK_DENS_;
    if (nwat_[a]>1) {
       TStrans_dw_[a] = -GASCONSTANT*BULK_DENS_*0.3*0.239*g_[a]*log(g_[a]);
       TStrans_norm_[a] = TStrans_dw_[a]/dens_[a];
       TStranstot_ += TStrans_dw_[a];
    }
    else {
       TStrans_dw_[a]=0; TStrans_norm_[a]=0;
    }
  }
  TStranstot_ *= Vvox_;
  mprintf("Total translational entropy of the grid: S_trans = %9.5f kcal/mol, Nf=%d\n", TStranstot_, NFRAME_);

  // Compute average voxel energy
  for (int a=0; a<MAX_GRID_PT_; a++) {
    if (nwat_[a]>1) {
       dEwh_dw_[a] = (wh_evdw_[a]+wh_eelec_[a])/(NFRAME_*Vvox_);
       //writeerror
       dEww_dw_ref_[a] = ((ww_evdw_[a]+ww_eelec_[a]) - nwat_[a]*BULK_E_)/(NFRAME_*Vvox_);
       dEwh_norm_[a] = (wh_evdw_[a]+wh_eelec_[a])/nwat_[a];
       //write error
       dEww_norm_ref_[a] = (ww_evdw_[a]+ww_eelec_[a])/nwat_[a] - BULK_E_;  
    }
    else {
       dEwh_dw_[a]=0; dEww_dw_ref_[a]=0; dEwh_norm_[a]=0; dEww_norm_ref_[a]=0;
    }
  }

  // Compute the actual voxel coordinates
  int voxel=0;
  for (int i = 0; i < griddim_[0]; ++i) {
    for (int j = 0; j < griddim_[1]; ++j) {
      for (int k = 0; k < griddim_[2]; ++k) {
	grid_x_[voxel] = Xcrd(i);
	grid_y_[voxel] = Xcrd(j);
	grid_z_[voxel] = Xcrd(k);
	voxel++;          
      }
    }
  }

  // Print the gist info file
  // Print the energy data
  if (!datafile_.empty()) {
    // Now write the data file with all of the GIST energies
    DataFile dfl;
    ArgList dummy;
    dfl.SetupDatafile(datafile_, dummy, 0);
    for (int i = 0; i < myDSL_.size(); i++) {
      dfl.AddSet(myDSL_[i]);
    }

    dfl.Write();

  }
  PrintDX("gist-dx");
  PrintOutput("gist-output");
}

  // Print GIST data in dx format
void Action_Gist::PrintDX(std::string const& filename)
{
  CpptrajFile outfile;
  if (outfile.OpenWrite(filename)) {
    mprinterr("Error: Could not open OpenDX output file.\n");
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
    outfile.Printf("%g %g %g\n", dens_[i], dens_[i+1], dens_[i+2]);
  // Print out any points we may have missed
  switch (MAX_GRID_PT_ % 3) {
    case 2: outfile.Printf("%g %g\n", dens_[MAX_GRID_PT_-2], dens_[MAX_GRID_PT_-1]); break;
    case 1: outfile.Printf("%g\n", dens_[MAX_GRID_PT_-1]); break;
  }

  outfile.CloseFile();
}

  // Print GIST data in dx format
void Action_Gist::PrintOutput(std::string const& filename)
{
  /*CpptrajFile outfile;
  if (outfile.OpenWrite(filename)) {
    mprinterr("Error: Could not open Gist output file.\n");
    return;
  }
  // Print the Output header
  outfile.Printf("GIST Output, information printed per voxel\n");
  outfile.Printf("xcoord ycoord zcoord population density normalized-density S-trans S-trans-normalized S-rot S-rot-normalized E-water-solute E-water-solute-normalized E-water-water E-water-water-normalized \n");
  
  // Now print out the data. 
  for (int i = 0; i < MAX_GRID_PT_ ; i++){
    outfile.Printf(" %10.5f %10.5f %10.5f\n %10.5f %10.5f %10.5f\n %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
		   grid_x_[i] , grid_y_[i], grid_z_[i], nwat_[i], dens_[i],g_[i], TStrans_dw_[i], TStrans_norm_[i], TSwNN_[i], TSNN_[i], dEwh_dw_[i],dEwh_norm_[i],dEww_dw_ref_[i],dEww_norm_ref_[i]);
  }
  
  outfile.CloseFile();
  */
  ofstream myfile;
  myfile.open("gist-output");  
  myfile <<"GIST Output, information printed per voxel\n";
  myfile <<"xcoord ycoord zcoord population density normalized-density S-trans S-trans-normalized S-rot S-rot-normalized E-water-solute E-water-solute-normalized E-water-water E-water-water-normalized \n";
  // Now print out the data. 
  for (int i = 0; i < MAX_GRID_PT_ ; i++){
    myfile << grid_x_[i] << " " << grid_y_[i]<< " " << grid_z_[i]<< " " << nwat_[i]<< " " << dens_[i]<< " " <<g_[i]<< " " << TStrans_dw_[i]<< " " << TStrans_norm_[i]<< " " << TSwNN_[i]<< " " << TSNN_[i]<< " " << dEwh_dw_[i]<< " " <<dEwh_norm_[i]<< " " <<dEww_dw_ref_[i]<< " " <<dEww_norm_ref_[i] << std::endl;
  }
  myfile.close();
}

// DESTRUCTOR
Action_Gist::~Action_Gist() {
  //fprintf(stderr,"Gist Destructor.\n");
  if (griddim_!=0) delete[] griddim_;
//  if (gridindex_!=0) delete[] gridindex_;
}
