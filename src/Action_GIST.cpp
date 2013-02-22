// GIST 
#include <cmath>
#include "Action_GIST.h"
#include "CpptrajStdio.h"
#include "Constants.h" // RADDEG
#include "Action_Pairwise.h"

// CONSTRUCTOR
Action_GIST::Action_GIST() :
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
  gridcntr_[0] = -1;
  gridcntr_[1] = -1;
  gridcntr_[2] = -1;
  
  griddim_[0] = -1;
  griddim_[1] = -1;
  griddim_[2] = -1;
  
  gridspacn_ = 0;
 } 

void Action_GIST::Help() {
  mprintf("gist <watermodel>[{tip3p|tip4p|tip4pew}] [gridcntr <xval> <yval> <zval>] [griddim <xval> <yval> <zval>] [gridspacn <spaceval>] [out <filename>] \n");
  mprintf("\tCalculate GIST between water molecules in selected site \n");
}

// Action_GIST::init()
Action::RetType Action_GIST::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
				  DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get keywords
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  DataSet::scalarType stype = DataSet::UNDEFINED;
  stype = DataSet::GIST;

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
    //griddim_[0] = 0.0;
    //griddim_[1] = 0.0;
    //griddim_[2] = 0.0;
    //mprintf("\tGIST grid dimension: %5.3f %5.3f %5.3f\n", griddim_[0],griddim_[1],griddim_[2]);
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

// Action_GIST::setup()
/** Set GIST up for this parmtop. Get masks etc.
  */
Action::RetType Action_GIST::Setup(Topology* currentParm, Topology** parmAddress) {
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
  return Action::OK;  
}


// Action_GIST::action()
Action::RetType Action_GIST::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {

  mprintf("GIST Action \n");
  //calculating energy
  atom_eelec_.assign(CurrentParm_->Natom(), 0);
  atom_evdw_.assign(CurrentParm_->Natom(), 0);

  // Set Masks
  std::string refmask = ":WAT";
  Mask1_.SetMaskString(refmask );
  refmask = ":WAT@O";
  Mask2_.SetMaskString(refmask );

  if (CurrentParm_->SetupIntegerMask( Mask1_ )) return Action::ERR;
  if (CurrentParm_->SetupIntegerMask( Mask2_ )) return Action::ERR;

 
  mprintf("GIST Action    : Atoms in mask1 [%s] %d \n",Mask1_.MaskString(),Mask1_.Nselected());
  mprintf("GIST Action    : Atoms in mask2 [%s] %d \n",Mask2_.MaskString(),Mask2_.Nselected());
  //  Action_Pairwise::NonbondEnergy( currentFrame, CurrentParm_, Mask1_);
  NonbondEnergy2( currentFrame, CurrentParm_, Mask1_ , Mask2_ );
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

void Action_GIST::NonbondEnergy2(Frame *frameIn, Topology *parmIn, AtomMask &maskIn , AtomMask &maskIn2) {
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
	atom_evdw_[atom2] += delta2;
	deltatest = delta2;
	// Cumulative eelec - divide between both atoms
	delta2 = e_elec * 0.5;
	atom_eelec_[atom1] += delta2;
	atom_eelec_[atom2] += delta2;
	//	mprintf("GIST Action NONBONDE atom1 %d atom2 %d eelec %f vdW %f \n",atom1,atom2, deltatest, delta2);

        // ----------------------------------------
      } // END pair not excluded
    } // END Inner loop
  } // END Outer loop

}

/*how can I call from action_pairwise (i can add a second declaration with the two masks, 
I can also try to access the vectors atom_eelec_ and atom_evdw_ to get the E values to 
everything else. If I can just call NonbondEnergy from Action_Pairwise and read atom_eelec_
 and atom_evdw_, I can remove these two functions  void Action_GIST::NonbondEnergy(Frame 
*frameIn, Topology *parmIn, AtomMask &maskIn , AtomMask &maskIn2) static void GetLJparam(
Topology const& top, double& A, double& B, int atom1, int atom2) */

