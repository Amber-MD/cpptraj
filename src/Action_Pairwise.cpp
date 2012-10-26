// Action_Pairwise
#include <cmath> //sqrt
#include "Action_Pairwise.h"
#include "CpptrajStdio.h"
#include "Trajout.h"
#include "Constants.h" // ELECTOAMBER

// CONSTRUCTOR
Action_Pairwise::Action_Pairwise() {
  //fprintf(stderr,"Pairwise Con\n");
  nb_calcType = NORMAL;
  RefParm=NULL;
  CurrentParm_ = 0;
  RefFrame=NULL;
  kes = 1.0;
  ELJ=0;
  Eelec=0;
  cut_eelec=1.0;
  cut_eelec1=-1.0;
  cut_evdw=1.0;
  cut_evdw1=-1.0;
  N_ref_interactions=0;
} 

void Action_Pairwise::Help() {
  mprintf("pairwise [<name>] [<mask>] [out <filename>] [cuteelec <cute>] [cutevdw <cutv>]\n");
  mprintf("         [ref <reffilename> | refindex <ref#>] [cutout <cutmol2name>]\n");
}

// DESTRUCTOR
Action_Pairwise::~Action_Pairwise() {
  //fprintf(stderr,"Pairwise Destructor.\n");
  Eout.CloseFile();
}

// Action_Pairwise::init()
Action::RetType Action_Pairwise::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get Keywords
  std::string dataout = actionArgs.GetStringKey("out");
  std::string eout = actionArgs.GetStringKey("eout");
  std::string referenceName = actionArgs.GetStringKey("ref");
  int refindex = actionArgs.getKeyInt("refindex",-1);
  cut_eelec = actionArgs.getKeyDouble("cuteelec",1.0);
  cut_eelec1 = -cut_eelec;
  cut_evdw = actionArgs.getKeyDouble("cutevdw",1.0);
  cut_evdw1 = -cut_evdw;
  cutout_ = actionArgs.GetStringKey("cutout");
  
  // Get Masks
  Mask0.SetMaskString( actionArgs.GetMaskNext() );
  std::string refmask = actionArgs.GetMaskNext();
  if (!refmask.empty())
    RefMask.SetMaskString(refmask);
  else
    RefMask.SetMaskString( Mask0.MaskString() );

  // Datasets
  std::string ds_name = actionArgs.GetStringNext();
  ds_vdw = DSL->AddSet(DataSet::DOUBLE, ds_name, "EVDW");
  ds_elec= DSL->AddSet(DataSet::DOUBLE, ds_name, "EELEC");
  // Add datasets to data file list
  DFL->AddSetToFile(dataout, ds_vdw);
  DFL->AddSetToFile(dataout, ds_elec);

  // Get reference structure
  if (!referenceName.empty() || refindex!=-1) {
    // Attempt to get reference index by name
    if (!referenceName.empty())
      refindex = FL->FindName( referenceName );
    // Get reference frame by index
    RefFrame = FL->GetFrame(refindex);
    if (RefFrame==NULL) {
      mprinterr("    Error: Pairwise::init: Could not get reference index %i\n",refindex);
      return Action::ERR;
    }
    // Set reference parm
    RefParm = FL->GetFrameParm(refindex);
    // Set up reference mask
    if ( RefParm->SetupIntegerMask(RefMask) ) return Action::ERR;
    if (RefMask.None()) {
      mprinterr("    Error: Pairwise::init: No atoms selected in reference mask.\n");
      return Action::ERR;
    }
    // Set up nonbonded params for reference
    if ( (N_ref_interactions=SetupNonbondParm( RefMask, RefParm )) == -1 ) return Action::ERR;
    // Calculate energy for reference
    nb_calcType = SET_REF;
    NonbondEnergy(RefFrame, RefParm, RefMask);
    mprintf("DEBUG:\tReference ELJ= %12.4lf  Eelec= %12.4lf\n",ELJ,Eelec);
    mprintf("DEBUG:\tSize of reference eelec array: %u\n",ref_nonbondEnergy.size());
    mprintf("DEBUG:\tSize of reference evdw array: %u\n",ref_nonbondEnergy.size());
    nb_calcType = COMPARE_REF;
  }

  // Output for individual atom energy | dEnergy
  if (!eout.empty()) {
    if (Eout.OpenWrite(eout)) {
      mprinterr("Error: Pairwise: Could not set up file %s for eout.\n",eout.c_str());
      return Action::ERR;
    }
  }

  // Action Info
  mprintf("    PAIRWISE: Atoms in mask [%s].\n",Mask0.MaskString());
  if (!eout.empty())
    mprintf("\tEnergy info for each atom will be written to %s\n",eout.c_str());
  if (RefFrame!=NULL) 
    mprintf("\tReference index %i, mask [%s]\n",refindex, RefMask.MaskString());
  mprintf("\tEelec absolute cutoff: %12.4lf\n",cut_eelec);
  mprintf("\tEvdw absolute cutoff: %12.4lf\n",cut_evdw);
  if (!cutout_.empty())
    mprintf("\tAtoms satisfying cutoff will be printed to %s.eX.mol2\n",
            cutout_.c_str());
  
  return Action::OK;
}

// Action_Pairwise::SetupNonbondParm()
/** Set up the exclusion list based on the given mask and parm.
  * \return the total number of interactions, -1 on error.
  */
int Action_Pairwise::SetupNonbondParm(AtomMask &maskIn, Topology *ParmIn) {
  // Charge is in units of electron charge, distance is in angstroms, so 
  // the electrostatic prefactor should be 332. However, since the charges
  // in Topology have presumably been converted from Amber charge units
  // create a new charged array multiplied by 18.2223. This makes calcs with 
  // Amber-converted charges more accurate at the cost of making non-Amber 
  // charges less accurate.
  atom_charge.clear();
  atom_charge.reserve( ParmIn->Natom() );
  for (Topology::atom_iterator atom = ParmIn->begin(); atom != ParmIn->end(); ++atom)
    atom_charge.push_back( (*atom).Charge() * ELECTOAMBER );
  // Check if LJ parameters present - need at least 2 atoms for it to matter.
  if (ParmIn->Natom() > 1 && (ParmIn->LJA().empty() || ParmIn->LJB().empty())) {
    mprinterr("Error: Pairwise::setup(): Parm does not have LJ information.\n");
    return -1;
  }

  // Determine the actual number of pairwise interactions that will be calcd.
  // This is ((N^2 - N) / 2) - SUM[ #excluded atoms]
  int N_interactions = ((maskIn.Nselected() * maskIn.Nselected()) - maskIn.Nselected()) / 2;
  for (AtomMask::const_iterator at = maskIn.begin(); at != maskIn.end(); ++at)
    N_interactions -= (*ParmIn)[ *at ].Nexcluded();

  // DEBUG - Print total number of interactions for this parm
  mprintf("\t%i interactions for this parm.\n",N_interactions);

  // DEBUG - Print exclusion list for each atom
  /*for (unsigned int atom = 0; atom < exclusionList.size(); atom++) {
    mprintf("\t%8u:",atom + 1);
    for (std::vector<int>::iterator eat = exclusionList[atom].begin();
                                    eat != exclusionList[atom].end();
                                    eat++)
    {
      mprintf(" %i",*eat + 1);
    }
    mprintf("\n");
  }*/
  return N_interactions;
}

// Action_Pairwise::setup()
/** Set up mask, allocate memory for exclusion list.
  */
Action::RetType Action_Pairwise::Setup(Topology* currentParm, Topology** parmAddress) {
  // Set up mask
  if ( currentParm->SetupIntegerMask( Mask0 ) ) return Action::ERR;
  if (Mask0.None()) {
    mprintf("    Error: Pairwise::setup: Mask has no atoms.\n");
    return Action::ERR;
  }

  // Set up exclusion list and determine total # interactions.
  int N_interactions = SetupNonbondParm(Mask0, currentParm);

  // If comparing to a reference frame for atom-by-atom comparison make sure
  // the number of interactions is the same in reference and parm.
  if (nb_calcType==COMPARE_REF) {
    if (N_interactions != N_ref_interactions) {
      mprinterr(
        "Error: Pairwise: # reference interactions (%i) != # interactions for this parm (%i)\n",
        N_ref_interactions,N_interactions
      );
      return Action::ERR;
    }
  }
  // Set up cumulative energy arrays
  atom_eelec.clear();
  atom_eelec.resize(currentParm->Natom(), 0);
  atom_evdw.clear();
  atom_evdw.resize(currentParm->Natom(), 0);
  // Print pairwise info for this parm
  Mask0.MaskInfo();
  CurrentParm_ = currentParm;      
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

// Action_Pairwise::NonbondEnergy()
/** Calculate non-bonded energy using the nonbondParm array. The total
  * LJ (vdw) energy is put in ELJ, and the total Coulomb (elec) energy
  * is put in Eelec. Depending on the value of nb_calcType, each pair
  * energy is either compared to a reference, distributed over both atoms
  * evenly in the cumulative array, or reference values are set. If comparing
  * to a reference structure, pairs for which the energy difference exceeds
  * the cutoffs are printed.
  */
void Action_Pairwise::NonbondEnergy(Frame *frameIn, Topology *parmIn, AtomMask &maskIn) {
  double delta2, Acoef, Bcoef;
  std::vector<NonbondEnergyType>::iterator refpair;
  NonbondEnergyType refE;

  ELJ = 0;
  Eelec = 0;
  refpair = ref_nonbondEnergy.begin();
  // Loop over all atom pairs and set information
  AtomMask::const_iterator mask_end = maskIn.end();
  AtomMask::const_iterator mask_end1 = maskIn.end();
  --mask_end1;
  // Outer loop
  for (AtomMask::const_iterator maskatom1 = maskIn.begin();
                                  maskatom1 != mask_end1; 
                                  maskatom1++)
  {
    // Set up coord index for this atom
    int coord1 = (*maskatom1) * 3;
    // Set up exclusion list for this atom
    Atom::excluded_iterator excluded_atom = (*parmIn)[*maskatom1].excludedbegin();
    // Inner loop
    AtomMask::const_iterator maskatom2 = maskatom1;
    ++maskatom2;
    for (; maskatom2 != mask_end; maskatom2++) {
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
        double r2=1/rij2;
        double r6=r2*r2*r2;
        double r12=r6*r6;
        double f12=Acoef*r12; // A/r^12
        double f6=Bcoef*r6;   // B/r^6
        double e_vdw=f12-f6;  // (A/r^12)-(B/r^6)
        ELJ += e_vdw;
        // LJ Force 
        //force=((12*f12)-(6*f6))*r2; // (12A/r^13)-(6B/r^7)
        //scalarmult(f,JI,F);
        // Coulomb energy 
        double qiqj = atom_charge[*maskatom1] * atom_charge[*maskatom2];
        double e_elec=kes * (qiqj/rij);
        Eelec += e_elec;
        // Coulomb Force
        //force=e_elec/rij; // kes*(qiqj/r)*(1/r)
        //scalarmult(f,JI,F);

        // ----------------------------------------
        int atom1 = *maskatom1;
        int atom2 = *maskatom2;
        // 1 - Comparison to reference, cumulative dEnergy on atoms
        if (nb_calcType == COMPARE_REF) {
          // dEvdw
          double delta_vdw = (*refpair).evdw - e_vdw;
          // dEelec
          double delta_eelec = (*refpair).eelec - e_elec;
          // Output
          if (Eout.IsOpen()) {
            if (delta_vdw > cut_evdw || delta_vdw < cut_evdw1) {
              Eout.Printf("\tAtom %6i@%4s-%6i@%4s dEvdw= %12.4lf\n",
                              atom1+1, (*parmIn)[atom1].c_str(),
                              atom2+1, (*parmIn)[atom2].c_str(), delta_vdw);
            }
            if (delta_eelec > cut_eelec || delta_eelec < cut_eelec1) {
              Eout.Printf("\tAtom %6i@%4s-%6i@%4s dEelec= %12.4lf\n",
                              atom1+1, (*parmIn)[atom1].c_str(),
                              atom2+1, (*parmIn)[atom2].c_str(),delta_eelec);
            }
          }
          // Divide the total pair dEvdw between both atoms.
          delta2 = delta_vdw * 0.5;
          atom_evdw[atom1] += delta2;
          atom_evdw[atom2] += delta2;
          // Divide the total pair dEelec between both atoms.
          delta2 = delta_eelec * 0.5;
          atom_eelec[atom1] += delta2;
          atom_eelec[atom2] += delta2;
        // 2 - No reference, just cumulative Energy on atoms
        } else if (nb_calcType == NORMAL) {
          // Cumulative evdw - divide between both atoms
          delta2 = e_vdw * 0.5;
          atom_evdw[atom1] += delta2;
          atom_evdw[atom2] += delta2;
          // Cumulative eelec - divide between both atoms
          delta2 = e_elec * 0.5;
          atom_eelec[atom1] += delta2;
          atom_eelec[atom2] += delta2;
        // 3 - Store the reference nonbond energy for this pair
        } else { // if nb_calcType == SET_REF
          refE.evdw = e_vdw;
          refE.eelec = e_elec;
          ref_nonbondEnergy.push_back( refE );
        }
        ++refpair;
        // ----------------------------------------
      } // END pair not excluded
    } // END Inner loop
  } // END Outer loop

}

// Action_Pairwise::WriteCutFrame()
int Action_Pairwise::WriteCutFrame(int frameNum, Topology *Parm, AtomMask& CutMask, 
                                   std::vector<double> const& CutCharges,
                                   Frame *frame, std::string const& outfilename) 
{
  Frame CutFrame(*frame, CutMask);
  // TEST: Write file containing only cut atoms
  Topology* CutParm = Parm->modifyStateByMask(CutMask);
  if (CutParm->Natom() != (int)CutCharges.size()) {
    mprinterr("Error: Pairwise: WriteCutFrame: # of charges (%u) != # mask atoms (%i)\n",
              CutCharges.size(), CutParm->Natom());
    delete CutParm;
    return 1;
  }
  std::vector<double>::const_iterator Qi = CutCharges.begin();
  for (Topology::iterator atom = Parm->begin(); atom != Parm->end(); ++atom)
    (*atom).SetCharge( *(Qi++) );
  //CutFrame.SetupFrame(CutParm->natom, CutParm->mass);
  //CutFrame.SetFrameFromMask(frame, CutMask);
  Trajout tout;
  if (tout.SetupTrajWriteWithArgs(outfilename,"multi",CutParm,TrajectoryFile::MOL2FILE)) {
    mprinterr("Error: Pairwise: Could not set up cut mol2 file %s\n",outfilename.c_str());
    delete CutParm;
    return 1;
  }
  tout.WriteFrame(frameNum,CutParm,CutFrame);
  tout.EndTraj();
  delete CutParm;
  return 0;
}

// Action_Pairwise::PrintCutAtoms()
/** Print atoms for which the cumulative energy satisfies the given
  * cutoffs. Also create MOL2 files containing those atoms.
  */
void Action_Pairwise::PrintCutAtoms(Frame *frame, int frameNum) {
  AtomMask CutMask; // TEST
  std::vector<double> CutCharges; // TEST
  // EVDW
  if (Eout.IsOpen()) {
    if (nb_calcType==COMPARE_REF)
      Eout.Printf("\tPAIRWISE: Cumulative dEvdw:");
    else
      Eout.Printf("\tPAIRWISE: Cumulative Evdw:");
    Eout.Printf(" Evdw < %.4lf, Evdw > %.4lf\n",cut_evdw1,cut_evdw);
  }
  for (int atom = 0; atom < CurrentParm_->Natom(); atom++) {
    if (atom_evdw[atom]>cut_evdw || atom_evdw[atom]<cut_evdw1) {
      if (Eout.IsOpen()) 
        Eout.Printf("\t\t%6i@%s: %12.4lf\n",atom+1,
                    (*CurrentParm_)[atom].c_str(),atom_evdw[atom]);
      CutMask.AddAtom(atom);
      CutCharges.push_back(atom_evdw[atom]);
    }
  }
  if (!cutout_.empty() && !CutMask.None()) {
    if (WriteCutFrame(frameNum, CurrentParm_, CutMask, CutCharges, frame, cutout_ + ".evdw.mol2")) 
      return;
  }
  CutMask.ResetMask();
  CutCharges.clear();
  // EELEC
  if (Eout.IsOpen()) {
    if (nb_calcType==COMPARE_REF)
      Eout.Printf("\tPAIRWISE: Cumulative dEelec:");
    else
      Eout.Printf("\tPAIRWISE: Cumulative Eelec:");
    Eout.Printf(" Eelec < %.4lf, Eelec > %.4lf\n",cut_eelec1,cut_eelec);
  }
  for (int atom = 0; atom < CurrentParm_->Natom(); atom++) { 
    if (atom_eelec[atom]>cut_eelec || atom_eelec[atom]<cut_eelec1) {
      if (Eout.IsOpen()) 
        Eout.Printf("\t\t%6i@%s: %12.4lf\n",atom+1,
                    (*CurrentParm_)[atom].c_str(), atom_eelec[atom]);
      CutMask.AddAtom(atom);
      CutCharges.push_back(atom_eelec[atom]);
    }  
  }
  if (!cutout_.empty() && !CutMask.None()) {
    if (WriteCutFrame(frameNum, CurrentParm_, CutMask, CutCharges, frame, cutout_ + ".eelec.mol2"))
      return;
  }
}

// Action_Pairwise::action()
Action::RetType Action_Pairwise::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) 
{
  //if (Energy(&Mask0, currentFrame, currentParm)) return 1;
  if (Eout.IsOpen()) Eout.Printf("PAIRWISE: Frame %i\n",frameNum);
  NonbondEnergy( currentFrame, CurrentParm_, Mask0 );
  PrintCutAtoms( currentFrame, frameNum );
  // Reset cumulative energy arrays
  atom_eelec.assign(CurrentParm_->Natom(), 0);
  atom_evdw.assign(CurrentParm_->Natom(), 0);

  ds_vdw->Add(frameNum, &ELJ);
  ds_elec->Add(frameNum, &Eelec);

  return Action::OK;
} 

// Action_Pairwise::print()
void Action_Pairwise::Print() {
/*  if (RefFrame!=NULL) {
    mprintf("\tPAIRWISE: Cumulative dEelec:\n");
    int iatom = 0;
    for (std::vector<double>::iterator atom = atom_eelec.begin();
         atom != atom_eelec.end(); atom++)
    {
      mprintf("\t\t%6i: %12.4lf\n",iatom++,*atom);
    }
  }*/
}
