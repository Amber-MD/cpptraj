// Gist 
#include <cmath>
#include <iostream> // cout
#include "Action_Gist.h"
#include "CpptrajFile.h"
#include "CpptrajStdio.h"
#include "Constants.h" // RADDEG
#include "DistRoutines.h"
#include "DataSet_integer.h"
#include "Box.h"
//#include "Topology.cpp"

// CONSTRUCTOR
Action_Gist::Action_Gist() :
  gist_(0),
  CurrentParm_(0),
  kes_(1.0),
  ELJ_(0),
  Eelec_(0),
  watermodel_(false),
  useTIP3P_(false),
  useTIP4P_(false),
  useTIP4PEW_(false)
{
  mprintf("\tGIST: INIT \n");
  gridcntr_[0] = -1;
  gridcntr_[1] = -1;
  gridcntr_[2] = -1;
  
  griddim_[0] = -1;
  griddim_[1] = -1;
  griddim_[2] = -1;
  
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
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  DataSet::scalarType stype = DataSet::UNDEFINED;
  stype = DataSet::GIST;

    mprintf("\tGIST: init2 \n");

  useTIP3P_ = actionArgs.hasKey("tip3p");
  useTIP4P_ = actionArgs.hasKey("tip4p");
  useTIP4PEW_ = actionArgs.hasKey("tip4pew");
  if (!useTIP3P_ && !useTIP4P_ && !useTIP4PEW_) {
    mprinterr("Error: gist: Only water models supprted are TIP3P and TIP4P\n");
    return Action::ERR;
  }
  
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

  
  if ( actionArgs.hasKey("griddim") ){
    griddim_[0] = actionArgs.getNextDouble(-1);
    griddim_[1] = actionArgs.getNextDouble(-1);
    griddim_[2] = actionArgs.getNextDouble(-1);
    mprintf("\tGIST grid dimension: %5.3f %5.3f %5.3f\n", griddim_[0],griddim_[1],griddim_[2]);
  }
  else{
    mprintf("\tGIST: No grid dimensiom values were found, using default (box size) \n");
    griddim_[0] = 30.0;
    griddim_[1] = 30.0;
    griddim_[2] = 30.0;
    mprintf("\tGIST grid dimension: %5.3f %5.3f %5.3f\n", griddim_[0],griddim_[1],griddim_[2]);
  }


  gridspacn_ = actionArgs.getKeyDouble("gridspacn", 0.50);
  mprintf("\tGIST grid spacing: %5.3f \n", gridspacn_);


  // Dataset to store gist results
  gist_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), "Gist");
  if (gist_==0) return Action::ERR;
  gist_->SetScalar( DataSet::M_DISTANCE, stype );
  // Add dataset to data file
  if (outfile != 0) outfile->AddSet( gist_ );
 
  return Action::OK;
}

// Action_Gist::setup()
/** Set Gist up for this parmtop. Get masks etc.
  */
Action::RetType Action_Gist::Setup(Topology* currentParm, Topology** parmAddress) {
  mprintf("GIST Setup \n");

  CurrentParm_ = currentParm;      
  // Set up cumulative energy arrays
  atom_eelec_.clear();
  atom_eelec_.resize(currentParm->Natom(), 0);
  atom_evdw_.clear();
  atom_evdw_.resize(currentParm->Natom(), 0);
  atom_charge_.clear();
  atom_charge_.reserve( currentParm->Natom() );
  for (Topology::atom_iterator atom = currentParm->begin(); atom != currentParm->end(); ++atom)
    atom_charge_.push_back( (*atom).Charge() * ELECTOAMBER );
  gridwat_.clear();
  gridwat_.reserve( currentParm->Natom() );

  // Set Masks
  std::string refmask = ":WAT";
  Mask1_.SetMaskString(refmask );
  refmask = ":WAT@O";
  Mask2_.SetMaskString(refmask );
  //refmask = ":!WAT";
  //Mask3_.SetMaskString(refmask );

  if (CurrentParm_->SetupIntegerMask( Mask1_ )) return Action::ERR;
  if (CurrentParm_->SetupIntegerMask( Mask2_ )) return Action::ERR;
  //if (CurrentParm_->SetupIntegerMask( Mask3_ )) return Action::ERR;

 
  mprintf("GIST Setup    : Atoms in mask1 [%s] %d \n",Mask1_.MaskString(),Mask1_.Nselected());
  mprintf("GIST Setup    : Atoms in mask2 [%s] %d \n",Mask2_.MaskString(),Mask2_.Nselected());
  //mprintf("GIST Setup    : Atoms in mask3 [%s] %d \n",Mask3_.MaskString(),Mask3_.Nselected());


  // Set up grid origin
  gridorig_[0] = gridcntr_[0] - 0.5*griddim_[0]*gridspacn_;
  gridorig_[1] = gridcntr_[1] - 0.5*griddim_[1]*gridspacn_;
  gridorig_[2] = gridcntr_[2] - 0.5*griddim_[2]*gridspacn_;
  mprintf("\tGIST grid origin: %5.3f %5.3f %5.3f\n", gridorig_[0],gridorig_[1],gridorig_[2]);
 
  MAX_GRID_PT=0.0;
  return Action::OK;  
}


// Action_Gist::action()
Action::RetType Action_Gist::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {

  mprintf("GIST Action \n");
  //calculating energy
  atom_eelec_.assign(CurrentParm_->Natom(), 0);
  atom_evdw_.assign(CurrentParm_->Natom(), 0);
  gridwat_.assign(CurrentParm_->Natom(), 0);

  //select water molecules
  Grid(currentFrame,  CurrentParm_);
  NonbondEnergy2( currentFrame, CurrentParm_, Mask1_ , Mask2_ );
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

void Action_Gist::NonbondEnergy2(Frame *frameIn, Topology *parmIn, AtomMask &maskIn , AtomMask &maskIn2) {
  double delta2, Acoef, Bcoef, deltatest;

  mprintf("GIST NonbondEnergy2  \n");
  mprintf("GIST    NonbondEnergy2: Atoms in mask1 [%s] %d \n",maskIn.MaskString(),maskIn.Nselected());
  mprintf("GIST    NonbondEnergy2: Atoms in mask2 [%s] %d \n",maskIn2.MaskString(),maskIn2.Nselected());

  ELJ_ = 0.0;
  Eelec_ = 0.0;
  // Loop over all atom pairs and set information
  AtomMask::const_iterator mask_end = maskIn.end();
  AtomMask::const_iterator mask_end2 = maskIn2.end();
  //  --mask_end;
  //--mask_end2;
  // Outer loop
  for (AtomMask::const_iterator maskatom1 = maskIn.begin();
                                  maskatom1 != mask_end; 
                                  maskatom1++)
  {
    // Set up coord index for this atom
    int coord1 = (*maskatom1) * 3;
    // Set up exclusion list for this atom
    Atom::excluded_iterator excluded_atom = (*parmIn)[*maskatom1].excludedbegin();
    // Inner loop
    //    AtomMask::const_iterator maskatom2 = maskatom1;
    
    //    ++maskatom2;
   for (AtomMask::const_iterator maskatom2 = maskIn2.begin(); 
	maskatom2 != mask_end2; 
	maskatom2++) {
      // If atom is excluded, just increment to next excluded atom;
      // otherwise perform energy calc.
      if ( excluded_atom != (*parmIn)[*maskatom1].excludedend() && *maskatom2 == *excluded_atom )
        ++excluded_atom;
      else {
        // Set up coord index for this atom
        int coord2 = (*maskatom2) * 3;
        // Calculate the vector pointing from atom2 to atom1
        Vec3 JI = Vec3(frameIn->CRD(coord1)) - Vec3(frameIn->CRD(coord2));
        double rij2 = JI.Magnitude2();
        // Normalize
        double rij = sqrt(rij2);
        JI /= rij;
        // LJ energy 
        GetLJparam(*parmIn, Acoef, Bcoef, *maskatom1, *maskatom2);
        double r2    = 1 / rij2;
        double r6    = r2 * r2 * r2;
        double r12   = r6 * r6;
        double f12   = Acoef * r12;  // A/r^12
        double f6    = Bcoef * r6;   // B/r^6
        double e_vdw = f12 - f6;     // (A/r^12)-(B/r^6)
        ELJ_ += e_vdw;
        // LJ Force 
        //force=((12*f12)-(6*f6))*r2; // (12A/r^13)-(6B/r^7)
        //scalarmult(f,JI,F);
        // Coulomb energy 
        double qiqj = atom_charge_[*maskatom1] * atom_charge_[*maskatom2];
        double e_elec = kes_ * (qiqj/rij);
        Eelec_ += e_elec;
        // Coulomb Force
        //force=e_elec/rij; // kes_*(qiqj/r)*(1/r)
        //scalarmult(f,JI,F);

        // ----------------------------------------
        int atom1 = *maskatom1;
        int atom2 = *maskatom2;
        
	// Cumulative evdw - divide between both atoms
	delta2 = e_vdw * 0.5;
	atom_evdw_[atom1] += delta2;
	//	atom_evdw_[atom2] += delta2;
	deltatest = delta2;
	// Cumulative eelec - divide between both atoms
	delta2 = e_elec * 0.5;
	atom_eelec_[atom1] += delta2;
	//	atom_eelec_[atom2] += delta2;
	//	mprintf("GIST Action NONBONDE atom1 %d atom2 %d eelec %f vdW %f \n",atom1,atom2, deltatest, delta2);

        // ----------------------------------------
      } // END pair not excluded
    } // END Inner loop
  } // END Outer loop

}


// Action_Gist::Grid()
void Action_Gist::Grid(Frame *frameIn, Topology* CurrentParm) {
  
  resnum=0;
  for (Topology::mol_iterator solvmol = CurrentParm_->MolStart();
       solvmol != CurrentParm_->MolEnd(); ++solvmol)
    {
      if (!(*solvmol).IsSolvent()) continue;
      int i = (*solvmol).BeginAtom();
      Vec3 O_wat = Vec3(frameIn->XYZ(i));
      gridwat_[resnum] = 10000000;
      // get the components of the water vector
      Vec3 comp = Vec3(O_wat) - Vec3(gridcntr_);
      double rij = sqrt( comp.Magnitude2() );
      Vec3 compnew = comp/gridspacn_;
      Vec3 index;
      index[0] = (int) compnew[0];
      index[1] = (int) compnew[1];
      index[2] = (int) compnew[2];
      if (index[0]>=0 && index[1]>=0 && index[2]>=0 && (index[0]<griddim_[0]) && (index[1]<griddim_[1]) && (index[2]<griddim_[2]))
	{
	  // this water belongs to grid point index[0], index[1], index[2]
	  int voxel = (index[0]*griddim_[1] + index[1])*griddim_[2] + index[2];
	  gridwat_[resnum] = voxel;
	}
      resnum++;
    }
  
}

void Action_Gist::EulerAngle(Frame *frameIn, Topology* CurrentParm) {

  Vec3 O_wat, H1_wat, H2_wat;
  
  //select water molecules
  int solventMolecules_ = CurrentParm_->Nsolvent();
  resnum=0;
  for (Topology::mol_iterator solvmol = CurrentParm_->MolStart();
       solvmol != CurrentParm_->MolEnd(); ++solvmol)
    {
      if (!(*solvmol).IsSolvent()) continue;

      if (gridwat_[resnum]>=MAX_GRID_PT) continue;

      int i = (*solvmol).BeginAtom();
      O_wat = Vec3(frameIn->XYZ(i));
      i++;
      H1_wat = Vec3(frameIn->XYZ(i));
      i++;
      H2_wat = Vec3(frameIn->XYZ(i));
      
      
      // Define lab frame of reference
      x_lab[0]=1.0; x_lab[1]=0; x_lab[2]=0;
      y_lab[0]=0; y_lab[1]=1.0; y_lab[2]=0;
      z_lab[0]=0; z_lab[1]=0; z_lab[2]=1.0;     
      
      // Define the water frame of reference - all axes must be normalized
      // make h1 the water x-axis (but first need to normalized)
      Vec3 x_wat = H1_wat;
      double rR_h1 = H1_wat.Normalize();
      // the normalized z-axis is the cross product of h1 and h2 
      Vec3 z_wat = x_wat.Cross( H2_wat );
      double rR_z_wat = z_wat.Normalize();
      // make y-axis as the cross product of h1 and z-axis
      Vec3 y_wat = z_wat.Cross( x_wat );
      double rR_y_wat = z_wat.Normalize();
      
      // Find the X-convention Z-X'-Z'' Euler angles between the water frame and the lab/host frame
      // First, theta = angle between the water z-axis of the two frames
      double dp = z_lab*( z_wat);
      double theta, phi, psi;
      theta = acos(dp);
      if (theta>0 && theta<PI) {
	// phi = angle between the projection of the water x-axis and the node
	// line of node is where the two xy planes meet = must be perpendicular to both z axes
	// direction of the lines of node = cross product of two normals (z axes)
	// acos of x always gives the angle between 0 and pi, which is okay for theta since theta ranges from 0 to pi
	Vec3 node = z_lab.Cross( z_wat );
	// Second, find the angle phi, which is between x_lab and the node
	dp = node*( x_lab );
        if (dp<=-1.0) phi = PI;
	else if (dp>=1.0) phi = PI;
	else phi = acos(dp);
	// check angle phi
	if (phi>0 && phi<(2*PI)) {
	  // method 2
	  Vec3 v = x_lab.Cross( node );
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
	/*Vec3 x_res = x_wat*( mat_W );
	Vec3 y_res = y_wat*( mat_W );
	Vec3 z_res = z_wat*( mat_W );*/
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
/*	rRx = x_res[0]*x_lab[0] + x_res[1]*x_lab[1] + x_res[2]*x_lab[2];
	rRy = y_res[0]*y_lab[0] + y_res[1]*y_lab[1] + y_res[2]*y_lab[2];
	rRz = z_res[0]*z_lab[0] + z_res[1]*z_lab[1] + z_res[2]*z_lab[2];
*/


/*	if (rRx>1+1E-6 || rRx<1-1E-6 || rRy>1+1E-6 || rRy<1-1E-6 || rRz>1+1E-6 || rRz<1-1E-6) {
	  std::cout << "fm=" << fm << ", wat=" << wat << ", gr=" << gr_pt << ", k=" << k << " ROTATION IS BAD!" << std::endl;
	  std::cout << "rx=" << rRx << ", ry=" << rRy << ", rz=" << rRz << std::endl;
	  std::cout << "water new x axis: " << x_res[0] << " " << x_res[1] << " " << x_res[2] << std::endl;
	  std::cout << "water new y axis: " << y_res[0] << " " << y_res[1] << " " << y_res[2] << std::endl;
	  std::cout << "water new z axis: " << z_res[0] << " " << z_res[1] << " " << z_res[2] << std::endl;
	  exit(EXIT_FAILURE);
	}	
	
	if (!(theta<=PI && theta>=0 && phi<=2*PI && phi>=0 && psi<=2*PI && psi>=0)) {
	  std::cout << "angles don't fall into range " << theta << " " << phi << " " << psi << std::endl;
	  std::cout << h1_wat[0] << " " << h1_wat[1] << " " << h1_wat[2] << " " << h2_wat[0] << " " << h2_wat[1] << " " << h2_wat[2] << std::endl;
	  exit(EXIT_FAILURE);
	}
*/
	int voxel = gridwat_[resnum];
	the_vox_[voxel].push_back(theta);
	phi_vox_[voxel].push_back(phi);
	psi_vox_[voxel].push_back(psi);
      }
      //     else std::cout << fm << " " << wat << " gimbal lock problem, two z_wat paralell" << std::endl;
      resnum++;
    }
  // DEBUG
  mprintf("GIST: Found %d solvent residues \n", resnum);
  if (solventMolecules_ != resnum) {
    mprinterr("GIST Error: No solvent molecules don't match %d %d\n", solventMolecules_, resnum);
  }
  
} 


