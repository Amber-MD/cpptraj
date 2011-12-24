// Pairwise
#include <cmath>
#include <cstdio>
#include "Action_Pairwise.h"
#include "CpptrajStdio.h"
#include "vectormath.h"
#include "TrajectoryFile.h"

// CONSTRUCTOR
Pairwise::Pairwise() {
  //fprintf(stderr,"Pairwise Con\n");
  nb_calcType = NORMAL;
  RefParm=NULL;
  RefFrame=NULL;
  kes = 1.0;
  ELJ=0;
  Eelec=0;
  cut_eelec=1.0;
  cut_eelec1=-1.0;
  cut_evdw=1.0;
  cut_evdw1=-1.0;
  cutout=NULL;
  N_ref_interactions=0;
} 

// DESTRUCTOR
Pairwise::~Pairwise() {
  //fprintf(stderr,"Pairwise Destructor.\n");
  Eout.CloseFile();
}

// Pairwise::init()
/** Expected call: pairwise [<name>] [<mask>] [out <filename>] [cuteelec <cute>] [cutevdw <cutv>]
                           [ref <reffilename> | refindex <ref#>] [cutout <cutmol2name>]
 */
// Dataset name will be the last arg checked for. Check order is:
//    1) Keywords
//    2) Masks
//    3) Dataset name
int Pairwise::init( ) {
  char *mask0, *dataout, *eout, *refmask, *referenceName, *ds_name;
  int refindex;

  // Get Keywords
  dataout = actionArgs.getKeyString("out",NULL);
  eout = actionArgs.getKeyString("eout",NULL);
  referenceName=actionArgs.getKeyString("ref",NULL);
  refindex=actionArgs.getKeyInt("refindex",-1);
  cut_eelec = actionArgs.getKeyDouble("cuteelec",1.0);
  cut_eelec1 = -cut_eelec;
  cut_evdw = actionArgs.getKeyDouble("cutevdw",1.0);
  cut_evdw1 = -cut_evdw;
  cutout = actionArgs.getKeyString("cutout",NULL);
  
  // Get Masks
  mask0 = actionArgs.getNextMask();
  //fprintf(stdout,"    Mask 1: %s\n",mask0);
  Mask0.SetMaskString(mask0);
  refmask = actionArgs.getNextMask();
  if (refmask!=NULL)
    RefMask.SetMaskString(refmask);
  else
    RefMask.SetMaskString(mask0);

  // Datasets
  ds_name = actionArgs.getNextString();
  ds_vdw = DSL->AddMulti(DOUBLE, ds_name, "EVDW");
  ds_elec= DSL->AddMulti(DOUBLE, ds_name, "EELEC");
  // Add datasets to data file list
  DFL->Add(dataout,ds_vdw);
  DFL->Add(dataout,ds_elec);

  // Get reference structure
  if (referenceName!=NULL || refindex!=-1) {
    // Attempt to get reference index by name
    if (referenceName!=NULL)
      refindex=FL->GetFrameIndex(referenceName);
    // Get reference frame by index
    RefFrame=FL->GetFrame(refindex);
    if (RefFrame==NULL) {
      mprinterr("    Error: Pairwise::init: Could not get reference index %i\n",refindex);
      return 1;
    }
    // Set reference parm
    RefParm=FL->GetFrameParm(refindex);
    // Set up reference mask
    if ( RefParm->SetupIntegerMask(RefMask, activeReference) ) return 1;
    if (RefMask.None()) {
      mprinterr("    Error: Pairwise::init: No atoms selected in reference mask.\n");
      return 1;
    }
    // Set up nonbonded params for reference
    if ( (N_ref_interactions=SetupNonbondParm( RefMask, RefParm )) == -1 ) return 1;
    // Calculate energy for reference
    nb_calcType = SET_REF;
    NonbondEnergy(RefFrame, RefParm, RefMask);
    mprintf("DEBUG:\tReference ELJ= %12.4lf  Eelec= %12.4lf\n",ELJ,Eelec);
    mprintf("DEBUG:\tSize of reference eelec array: %u\n",ref_nonbondEnergy.size());
    mprintf("DEBUG:\tSize of reference evdw array: %u\n",ref_nonbondEnergy.size());
    nb_calcType = COMPARE_REF;
  }

  // Output for individual atom energy | dEnergy
  if (eout!=NULL) {
    if (Eout.SetupFile(eout,WRITE,debug)) {
      mprinterr("Error: Pairwise: Could not set up file %s for eout.\n",eout);
      return 1;
    }
    Eout.OpenFile();
  }

  // Action Info
  mprintf("    PAIRWISE: Atoms in mask [%s].\n",Mask0.MaskString());
  if (eout!=NULL)
    mprintf("              Energy info for each atom will be written to %s\n",eout);
  if (RefFrame!=NULL) 
    mprintf("              Reference index %i, mask [%s]\n",refindex, RefMask.MaskString());
  mprintf("              Eelec absolute cutoff: %12.4lf\n",cut_eelec);
  mprintf("              Evdw absolute cutoff: %12.4lf\n",cut_evdw);
  if (cutout!=NULL)
    mprintf("              Atoms satisfying cutoff will be printed to %s.eX.mol2\n",cutout);
  
  return 0;
}

// Pairwise::SetupNonbondParm()
/** Set up the exclusion list based on the given mask and parm.
  * \return the total number of interactions, -1 on error.
  */
int Pairwise::SetupNonbondParm(AtomMask &maskIn, AmberParm *ParmIn) {
  int N_interactions = 0;

  // Charge is in units of electron charge, distance is in angstroms, so 
  // the electrostatic prefactor should be 332. However, since the charges
  // in AmberParm have presumably been converted from Amber charge units
  // create a new charged array multiplied by 18.2223. This makes calcs with 
  // Amber-converted charges more accurate at the cost of making non-Amber 
  // charges less accurate.
  if ( ParmIn->AmberChargeArray( atom_charge ) ) {
    mprinterr("Error: Pairwise::setup(): Parm does not have charge information.\n");
    return -1;
  }
  // Check if LJ parameters present - need at least 2 atoms for it to matter.
  if (ParmIn->natom>1) {
    double Atemp = 0;
    double Btemp = 0;
    if (ParmIn->GetLJparam(&Atemp, &Btemp, 0, 1)) {
      mprinterr("Error: Pairwise::setup(): Parm does not have LJ information.\n");
      return -1;
    }
  }

  // Set up exclusion list for each atom in mask. 
  N_interactions = ParmIn->SetupExcludedAtomsList(maskIn, exclusionList);

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

// Pairwise::setup()
/** Set up mask, allocate memory for exclusion list.
  */
int Pairwise::setup() {
  // Set up mask
  if ( currentParm->SetupIntegerMask( Mask0, activeReference) ) return 1;
  if (Mask0.None()) {
    mprintf("    Error: Pairwise::setup: Mask has no atoms.\n");
    return 1;
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
      return 1;
    }
  }
  // Set up cumulative energy arrays
  atom_eelec.clear();
  atom_eelec.resize(currentParm->natom, 0);
  atom_evdw.clear();
  atom_evdw.resize(currentParm->natom, 0);
  // Print pairwise info for this parm
  mprintf("    PAIRWISE: Mask %s corresponds to %i atoms.\n",Mask0.MaskString(), Mask0.Nselected);
        
  return 0;  
}

// Pairwise::NonbondEnergy()
/** Calculate non-bonded energy using the nonbondParm array. The total
  * LJ (vdw) energy is put in ELJ, and the total Coulomb (elec) energy
  * is put in Eelec. Depending on the value of nb_calcType, each pair
  * energy is either compared to a reference, distributed over both atoms
  * evenly in the cumulative array, or reference values are set. If comparing
  * to a reference structure, pairs for which the energy difference exceeds
  * the cutoffs are printed.
  */
void Pairwise::NonbondEnergy(Frame *frameIn, AmberParm *parmIn, AtomMask &maskIn) {
  double rij, rij2, JI[3], delta2, Acoef, Bcoef;
  std::vector<NonbondEnergyType>::iterator refpair;
  NonbondEnergyType refE;

  ELJ = 0;
  Eelec = 0;
  JI[0]=0; JI[1]=0; JI[2]=0;
  refpair = ref_nonbondEnergy.begin();
  // Loop over all atom pairs and set information
  std::vector<int>::iterator mask_end = maskIn.Selected.end();
  std::vector<int>::iterator mask_end1 = maskIn.Selected.end();
  --mask_end1;
  // Outer loop
  for (std::vector<int>::iterator maskatom1 = maskIn.Selected.begin();
                                  maskatom1 != mask_end1; 
                                  maskatom1++)
  {
    // Set up coord index for this atom
    int coord1 = (*maskatom1) * 3;
    // Set up exclusion list for this atom
    std::vector<int>::iterator excluded_atom = exclusionList[*maskatom1].begin();
    // Inner loop
    std::vector<int>::iterator maskatom2 = maskatom1;
    ++maskatom2;
    for (; maskatom2 != mask_end; maskatom2++) {
      // If atom is excluded, just increment to next excluded atom;
      // otherwise perform energy calc.
      if ( *maskatom2 == *excluded_atom )
        ++excluded_atom;
      else {
        // Set up coord index for this atom
        int coord2 = (*maskatom2) * 3;
        // Calculate the vector pointing from atom2 to atom1
        vector_sub(JI, frameIn->X+coord1, frameIn->X+coord2);
        // Normalize
        rij = vector_norm(JI, &rij2);
        // LJ energy 
        parmIn->GetLJparam(&Acoef, &Bcoef,*maskatom1,*maskatom2);
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
              Eout.IO->Printf("\tAtom %6i@%4s-%6i@%4s dEvdw= %12.4lf\n",
                              atom1+1,currentParm->AtomName(atom1),
                              atom2+1,currentParm->AtomName(atom2),delta_vdw);
            }
            if (delta_eelec > cut_eelec || delta_eelec < cut_eelec1) {
              Eout.IO->Printf("\tAtom %6i@%4s-%6i@%4s dEelec= %12.4lf\n",
                              atom1+1,currentParm->AtomName(atom1),
                              atom2+1,currentParm->AtomName(atom2),delta_eelec);
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

// Pairwise::WriteCutFrame()
int Pairwise::WriteCutFrame(AmberParm *Parm, AtomMask *CutMask, double *CutCharges,
                            Frame *frame, char *outfilename) 
{
  AmberParm *CutParm;
  Frame CutFrame;
  TrajectoryFile tout;
  // TEST: Write file containing only cut atoms
  CutParm = Parm->modifyStateByMask(CutMask->Selected, NULL);
  CutParm->SetCharges(CutCharges);
  CutFrame.SetupFrame(CutParm->natom, CutParm->mass);
  CutFrame.SetFrameFromMask(frame, CutMask);
  if (tout.SetupWriteWithArgs(outfilename,"multi",CutParm,MOL2FILE)) {
    mprinterr("Error: Pairwise: Could not set up cut mol2 file %s\n",outfilename);
    return 1;
  }
  tout.WriteFrame(frameNum,CutParm,CutFrame);
  tout.EndTraj();
  delete CutParm;
  return 0;
}

// Pairwise::PrintCutAtoms()
/** Print atoms for which the cumulative energy satisfies the given
  * cutoffs. Also create MOL2 files containing those atoms.
  */
void Pairwise::PrintCutAtoms(Frame *frame) {
  char buffer[256]; // NOTE: Temporary
  AtomMask CutMask; // TEST
  std::vector<double> CutCharges; // TEST
  // EVDW
  if (Eout.IsOpen()) {
    if (nb_calcType==COMPARE_REF)
      Eout.IO->Printf("\tPAIRWISE: Cumulative dEvdw:");
    else
      Eout.IO->Printf("\tPAIRWISE: Cumulative Evdw:");
    Eout.IO->Printf(" Evdw < %.4lf, Evdw > %.4lf\n",cut_evdw1,cut_evdw);
  }
  for (int atom = 0; atom < currentParm->natom; atom++) {
    if (atom_evdw[atom]>cut_evdw || atom_evdw[atom]<cut_evdw1) {
      if (Eout.IsOpen()) 
        Eout.IO->Printf("\t\t%6i@%s: %12.4lf\n",atom+1,currentParm->AtomName(atom),atom_evdw[atom]);
      CutMask.AddAtom(atom);
      CutCharges.push_back(atom_evdw[atom]);
    }
  }
  if (cutout!=NULL && !CutMask.None()) {
    sprintf(buffer,"%s.evdw.mol2",cutout);
    if (WriteCutFrame(currentParm, &CutMask, &CutCharges[0], frame, buffer)) return;
  }
  CutMask.ResetMask();
  CutCharges.clear();
  // EELEC
  if (Eout.IsOpen()) {
    if (nb_calcType==COMPARE_REF)
      Eout.IO->Printf("\tPAIRWISE: Cumulative dEelec:");
    else
      Eout.IO->Printf("\tPAIRWISE: Cumulative Eelec:");
    Eout.IO->Printf(" Eelec < %.4lf, Eelec > %.4lf\n",cut_eelec1,cut_eelec);
  }
  for (int atom = 0; atom < currentParm->natom; atom++) { 
    if (atom_eelec[atom]>cut_eelec || atom_eelec[atom]<cut_eelec1) {
      if (Eout.IsOpen()) 
        Eout.IO->Printf("\t\t%6i@%s: %12.4lf\n",atom+1,
                        currentParm->AtomName(atom),atom_eelec[atom]);
      CutMask.AddAtom(atom);
      CutCharges.push_back(atom_eelec[atom]);
    }  
  }
  if (cutout!=NULL && !CutMask.None()) {
    sprintf(buffer,"%s.eelec.mol2",cutout);
    if (WriteCutFrame(currentParm, &CutMask, &CutCharges[0], frame,buffer)) return;
  }
}

// Pairwise::action()
int Pairwise::action() {
  //if (Energy(&Mask0, currentFrame, currentParm)) return 1;
  if (Eout.IsOpen()) Eout.IO->Printf("PAIRWISE: Frame %i\n",frameNum);
  NonbondEnergy( currentFrame, currentParm, Mask0 );
  PrintCutAtoms( currentFrame );
  // Reset cumulative energy arrays
  atom_eelec.assign(currentParm->natom, 0);
  atom_evdw.assign(currentParm->natom, 0);

  ds_vdw->Add(frameNum, &ELJ);
  ds_elec->Add(frameNum, &Eelec);

  return 0;
} 

// Pairwise::print()
void Pairwise::print() {
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
