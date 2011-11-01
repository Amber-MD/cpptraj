// Pairwise
#include <cmath>
#include <cstdio>
#include "Action_Pairwise.h"
#include "CpptrajStdio.h"
#include "vectormath.h"
// TEST
//#include "Traj_Mol2File.h"
#include "TrajectoryFile.h"

// CONSTRUCTOR
Pairwise::Pairwise() {
  //fprintf(stderr,"Pairwise Con\n");
  RefParm=NULL;
  RefFrame=NULL;
  skipv = NULL;
  natexidx = NULL;
  hasExclusion=true;
  kes = 1.0;
  ELJ=0;
  Eelec=0;
  cut_eelec=1.0;
  cut_eelec1=-1.0;
  cut_evdw=1.0;
  cut_evdw1=-1.0;
  cutout=NULL;
} 

// DESTRUCTOR
Pairwise::~Pairwise() {
  //fprintf(stderr,"Pairwise Destructor.\n");
  if (skipv!=NULL) delete[] skipv;
  if (natexidx!=NULL) delete[] natexidx;
  Eout.CloseFile();
}

/* Pairwise::init()
 * Expected call: pairwise [<name>] [<mask>] [out <filename>] [cuteelec <cute>] [cutevdw <cutv>]
                           [ref <reffilename> | refindex <ref#>] [cutout <cutmol2name>]
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
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
    if (RefMask.SetupMask(RefParm, activeReference, debug)) return 1;
    if (RefMask.None()) {
      mprinterr("    Error: Pairwise::init: No atoms selected in reference mask.\n");
      return 1;
    }
    // Allocate exclusion list for reference
    AllocateExclusion(RefParm);
    // Set up comparison array
    N_ref_interactions = NumInteractions(&RefMask, RefParm);
    ref_eelec.clear();
    ref_eelec.resize( N_ref_interactions, 0 );
    ref_evdw.clear();
    ref_evdw.resize( N_ref_interactions, 0 );
    // Calculate energy for reference
    RefEnergy();
    mprintf("DEBUG:\tReference ELJ= %12.4lf  Eelec= %12.4lf\n",ELJ,Eelec);
    mprintf("DEBUG:\tSize of reference eelec array: %u\n",ref_eelec.size());
    mprintf("DEBUG:\tSize of reference evdw array: %u\n",ref_evdw.size());
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
  mprintf("    PAIRWISE: Atoms in mask [%s].\n",Mask0.maskString);
  if (eout!=NULL)
    mprintf("              Energy info for each atom will be written to %s\n",eout);
  if (RefFrame!=NULL) 
    mprintf("              Reference index %i, mask [%s]\n",refindex, RefMask.maskString);
  mprintf("              Eelec absolute cutoff: %12.4lf\n",cut_eelec);
  mprintf("              Evdw absolute cutoff: %12.4lf\n",cut_evdw);
  if (cutout!=NULL)
    mprintf("              Atoms satisfying cutoff will be printed to %s.eX.mol2\n",cutout);
  
  return 0;
}

/* Pairwise::AllocateExclusion()
 * Allocate memory for the exclusion list based on the given parm.
 */
int Pairwise::AllocateExclusion(AmberParm *Parm) {
  int idx = 0;

  if (skipv!=NULL) delete[] skipv;
  skipv = new bool[ Parm->natom ];
  for (int atom=0; atom < Parm->natom; atom++) 
    skipv[atom]=false;
  // Check if exclusion info present
  if (Parm->NumExcludedAtoms(0)==-1) hasExclusion=false;
  if (hasExclusion) {
    // Create an array holding indices for each atom into NATEX
    if (natexidx!=NULL) delete[] natexidx;
    natexidx = new int[ Parm->natom ];
    for (int atom=0; atom < Parm->natom; atom++) {
      natexidx[atom] = idx;
      idx += Parm->NumExcludedAtoms(atom);
    }
    // DEBUG
    if (debug>0) {
      mprintf("DEBUG: EXCLUDED ATOM IDX:\n");
      for (int atom=0; atom < Parm->natom; atom++) 
        mprintf("\tAtom %i: Index %i\n",atom,natexidx[atom]);
    }
  }
  // Check if LJ parameters present - need at least 2 atoms for it to matter.
  if (Parm->natom>1) {
    double Atemp = 0;
    double Btemp = 0;
    if (Parm->GetLJparam(&Atemp, &Btemp, 0, 1)) {
      mprinterr("Error: Pairwise::setup(): Parm does not have LJ information.\n");
      return 1;
    }
  }

  return 0;
}

/* Pairwise::NumInteractions()
 * Based on the given atom mask and parm determine the total number of
 * pairwise interactions that will be calculated.
 * There should be ((N^2 - N) / 2) - SUM(Numex) interactions, however
 * since the excluded atoms list can include 0 (which translates to 
 * atom -1, nonexistent) there may be more interactions than expected,
 * so each one must be explicitly counted.
 */
int Pairwise::NumInteractions(AtomMask *atommask, AmberParm *Parm) {
  int Nselected = atommask->Nselected;
  int N_interactions = 0;

  if (hasExclusion) {
    for (int maskidx = 0; maskidx < Nselected - 1; maskidx++) {
      int atom1 = atommask->Selected[maskidx];
      SetupExclusion(Parm, atom1);
      for (int maskidx2 = maskidx + 1; maskidx2 < Nselected; maskidx2++) {
        int atom2 = atommask->Selected[maskidx2];
        if (!skipv[atom2]) N_interactions++; 
      }
    }
  } else {
    N_interactions = (((Nselected * Nselected) - Nselected) / 2);
  }

  // DEBUG
  mprintf("DEBUG:\t%i total interactions.\n",N_interactions);
  return N_interactions;
}

/* Pairwise::setup()
 * Set up mask, allocate memory for exclusion list.
 */
int Pairwise::setup() {
  // Set up mask
  if ( Mask0.SetupMask(P,activeReference,debug) ) return 1;
  if (Mask0.None()) {
    mprintf("    Error: Pairwise::setup: Mask has no atoms.\n");
    return 1;
  }
  // Allocate memory for exclusion list
  if (AllocateExclusion(P)) return 1;
  // If a reference frame is defined for atom-by-atom comparison make sure 
  // the number of interactions is the same in reference and parm.
  if (RefFrame!=NULL) {
    int N_interactions = NumInteractions(&Mask0, P);
    if (N_interactions != N_ref_interactions) {
      mprinterr(
        "Error: Pairwise: # reference interactions (%i) != # interactions for this parm (%i)\n",
        N_ref_interactions, N_interactions
      );
      return 1;
    }
  }
  // Set up cumulative energy arrays
  atom_eelec.clear();
  atom_eelec.resize(P->natom, 0);
  atom_evdw.clear();
  atom_evdw.resize(P->natom, 0);
  // Print pairwise info for this parm
  mprintf("    PAIRWISE: Mask %s corresponds to %i atoms.\n",Mask0.maskString, Mask0.Nselected);
        
  return 0;  
}

/* Pairwise::SetupExclusion()
 * Sets up the exclusion list, skipv, for the given atom; skipv will be
 * true for atoms that will be skipped (i.e. excluded) in the non-bonded
 * calc. skipv should be allocd and natexidx should be set up prior to 
 * calling this routine.
 */
int Pairwise::SetupExclusion(AmberParm *Parm, int atom1) {
  int jexcl, jexcl_last, natex;
  // Set up exclusion list for atom1
  for (int atom2=atom1+1; atom2 < Parm->natom; atom2++) 
    skipv[atom2]=false;
  // Set start and stop indices into the Excluded atoms list, NATEX
  jexcl = natexidx[atom1];
  jexcl_last = jexcl + Parm->NumExcludedAtoms(atom1);
  for (int idx = jexcl; idx < jexcl_last; idx++) {
    natex = Parm->Natex(idx);
    if (debug>0) mprintf("\t\tPair %i - %i will be excluded.\n",atom1,natex);
    //if (natex >= P->natom) {
    //  mprinterr("Error: excluded atom index %i for atom %i > # atoms.\n",idx,atom1);
    //  return 1;
    //}
    if (natex>=0) skipv[ natex ] = true;
  }
  return 0;
}

/* Pairwise::Energy_LJ()
 * Calculate the Lennard-Jones 6-12 energy and force between two atoms
 * separated by the given distance squared.
 */
double Pairwise::Energy_LJ(AmberParm *Parm, int atom1, int atom2, double rij2, double *force) {
  double Acoef, Bcoef, r2, r6, r12, f12, f6, energy;
  // LJ Energy
  Parm->GetLJparam(&Acoef, &Bcoef,atom1,atom2);
  r2=1/rij2;
  r6=r2*r2*r2;
  r12=r6*r6;
  f12=Acoef*r12;               // A/r^12
  f6=Bcoef*r6;                 // B/r^6
  energy=f12-f6;               // (A/r^12)-(B/r^6)
  // LJ Force 
  *force=((12*f12)-(6*f6))*r2; // (12A/r^13)-(6B/r^7)
  //scalarmult(f,JI,F);
  //M[i].f[0]+=f[0];
  //M[i].f[1]+=f[1];
  //M[i].f[2]+=f[2];
  //M[j].f[0]-=f[0];
  //M[j].f[1]-=f[1];
  //M[j].f[2]-=f[2];
  if (debug>0) {
    mprintf("\t\t\tLJ Force:  A= %lf  B= %lf  F= %lf  E= %lf\n",Acoef,Bcoef,force,energy);
    //mprintf("\t\t\tevdw= %lf\n",ELJ);
    mprintf("\t\t\t1/r6= %lf f6=%lf 1/r12= %lf f12=%lf\n",r6,f6,r12,f12);
    //mprintf("\t\t\tAtom %i : Fx= %lf  Fy= %lf  Fz= %lf\n",atom1,f[0],f[1],f[2]);
    //mprintf("\t\t\tAtom %i : Fx= %lf  Fy= %lf  Fz= %lf\n",atom2,-f[0],-f[1],-f[2]);
  }
  return energy;
}

/* Pairwise::Energy_Coulomb()
 * Calculate the Coulomb electrostatic energy and force between two atoms 
 * separated by the given distance.
 */
double Pairwise::Energy_Coulomb(AmberParm *Parm, int atom1, int atom2, double rij, double *force) {
  double qi, qj, qiqj, energy;
  // Coulomb Energy 
  // NOTE: Currently cpptraj converts charge to units of e- when topology
  // is read in. Need to convert back here. Not really efficient, just a 
  // hack for now.
  qi = Parm->charge[atom1] * 18.2223;
  qj = Parm->charge[atom2] * 18.2223;
  qiqj = qi * qj;
  energy=kes * (qiqj/rij);
  // Coulomb Force
  *force=energy/rij; // kes*(qiqj/r)*(1/r)
  //scalarmult(f,JI,F);
  //M[i].f[0]+=f[0];
  //M[i].f[1]+=f[1];
  //M[i].f[2]+=f[2];
  //M[j].f[0]-=f[0];
  //M[j].f[1]-=f[1];
  //M[j].f[2]-=f[2];
  if (debug>0) {
    mprintf("\t\t\tCoulomb Force:  q%i= %lf  q%i= %lf  F= %lf  E= %lf\n",
            atom1,Parm->charge[atom1],atom2,Parm->charge[atom2],force,energy);
    //mprintf("\t\t\tAtom %i : Fx= %lf  Fy= %lf  Fz= %lf\n",atom1,f[0],f[1],f[2]);
    //mprintf("\t\t\tAtom %i : Fx= %lf  Fy= %lf  Fz= %lf\n",atom2,-f[0],-f[1],-f[2]);
  }
  return energy;
}

/* Pairwise::WriteCutFrame()
 */
int Pairwise::WriteCutFrame(AmberParm *Parm, AtomMask *CutMask, double *CutCharges,
                            Frame *frame, char *outfilename) 
{
  AmberParm *CutParm;
  Frame CutFrame;
  TrajectoryFile tout;
  // TEST: Write file containing only cut atoms
  CutParm = Parm->modifyStateByMask(CutMask->Selected, CutMask->Nselected);
  CutParm->SetCharges(CutCharges);
  CutFrame.SetupFrame(CutParm->natom, CutParm->mass);
  CutFrame.SetFrameFromMask(frame, CutMask);
  if (tout.SetupWriteWithArgs(outfilename,"multi",CutParm,MOL2FILE)) {
    mprinterr("Error: Pairwise: Could not set up cut mol2 file %s\n",outfilename);
    return 1;
  }
  tout.WriteFrame(currentFrame,CutParm,CutFrame.X, NULL, NULL, 0.0);
  tout.EndTraj();
  delete CutParm;
  return 0;
}

/* Pairwise::RefEnergy()
 * Calculate pairwise energy for the reference mask, frame, and parm.
 * Fill the ref_X arrays.
 */
void Pairwise::RefEnergy() {
  int atom1, atom2, Ncomparison;
  double rij, rij2, JI[3], force, e_vdw, e_elec;
  // Loop over all atom pairs excluding self
  ELJ = 0;
  Eelec = 0;
  JI[0]=0; JI[1]=0; JI[2]=0;
  Ncomparison = 0; // Pairwise interaction counter
  // Outer loop
  for (int maskidx1 = 0; maskidx1 < RefMask.Nselected - 1; maskidx1++) {
    atom1 = RefMask.Selected[maskidx1];
    if (hasExclusion)
      SetupExclusion(RefParm, atom1);
    // Inner loop
    for (int maskidx2 = maskidx1 + 1; maskidx2 < RefMask.Nselected; maskidx2++) {
      atom2 = RefMask.Selected[maskidx2];
      int coord1 = atom1 * 3;
      int coord2 = atom2 * 3;
      // Calculate the vector pointing from atom2 to atom1
      vector_sub(JI, RefFrame->X+coord1, RefFrame->X+coord2);
      // Normalize
      rij = vector_norm(JI, &rij2);
      // Non-bonded energy
      if (!skipv[atom2]) {
        // Lennard Jones 6-12 Energy
        e_vdw = Energy_LJ(RefParm, atom1, atom2, rij2, &force);
        ELJ += e_vdw;
        // Coulomb Energy
        e_elec = Energy_Coulomb(RefParm, atom1, atom2, rij, &force);
        Eelec += e_elec;
        // Store in ref_X arrays
        ref_evdw[ Ncomparison ] = e_vdw;
        ref_eelec[ Ncomparison ] = e_elec;
        Ncomparison++;
      }
    } // End inner loop
  } // End outer loop
}
      
/* Pairwise::Energy()
 * Calculate pairwise energy for the given mask, frame, and parm. Sets
 * ELJ and Eelec. 
 */
int Pairwise::Energy(AtomMask *atommask, Frame *frame, AmberParm *Parm) {
  int atom1, atom2, Ncomparison;
  double rij, rij2, JI[3], force, e_vdw, e_elec, delta, delta2;
  char buffer[256]; // NOTE: Temporary
  // Loop over all atom pairs excluding self
  ELJ = 0;
  Eelec = 0;
  JI[0]=0; JI[1]=0; JI[2]=0;
  Ncomparison = 0; // Pairwise interaction counter
  // Data output
  if (Eout.IsOpen())
    Eout.IO->Printf("PAIRWISE: Frame %i\n",currentFrame);
  // Outer loop
  for (int maskidx1 = 0; maskidx1 < atommask->Nselected - 1; maskidx1++) {
    atom1 = atommask->Selected[maskidx1];
    if (debug>0) mprintf("\tPAIRWISE: Atom %i\n",atom1+1);
    // Set up exclusion list if necessary
    if (hasExclusion) 
      SetupExclusion(Parm, atom1);    
    // Inner loop
    for (int maskidx2 = maskidx1 + 1; maskidx2 < atommask->Nselected; maskidx2++) {
      atom2 = atommask->Selected[maskidx2];
      if (debug>0) mprintf("\t\tCalc for Pair %i - %i:\n",atom1,atom2);
      int coord1 = atom1 * 3;
      int coord2 = atom2 * 3;
      // Calculate the vector pointing from atom2 to atom1
      vector_sub(JI, frame->X+coord1, frame->X+coord2);
      // Normalize
      rij = vector_norm(JI, &rij2);
      // Non-bonded energy
      if (!skipv[atom2]) {
        // Lennard Jones 6-12 Energy
        e_vdw = Energy_LJ(Parm, atom1, atom2, rij2, &force);
        ELJ += e_vdw;
        // Coulomb Energy
        e_elec = Energy_Coulomb(Parm, atom1, atom2, rij, &force);
        Eelec += e_elec;
        // 1 - Comparison to reference, cumulative dEnergy on atoms
        if (RefFrame!=NULL) {
          // dEvdw
          delta = ref_evdw[ Ncomparison ] - e_vdw;
          if (Eout.IsOpen() && (delta > cut_evdw || delta < cut_evdw1)) {
            Eout.IO->Printf("\tAtom %6i@%4s-%6i@%4s dEvdw= %12.4lf\n",atom1+1,Parm->names[atom1],
                            atom2+1,Parm->names[atom2],delta);
          }
          // Divide the total pair dEvdw between both atoms.
          delta2 = delta * 0.5;
          atom_evdw[atom1] += delta2;
          atom_evdw[atom2] += delta2;
          // dEelec
          delta = ref_eelec[ Ncomparison ] - e_elec;
          if (Eout.IsOpen() && (delta > cut_eelec || delta < cut_eelec1)) {
            Eout.IO->Printf("\tAtom %6i@%4s-%6i@%4s dEelec= %12.4lf\n",atom1+1,Parm->names[atom1],
                            atom2+1,Parm->names[atom2],delta);
          }
          // Divide the total pair dEelec between both atoms.
          delta2 = delta * 0.5;
          atom_eelec[atom1] += delta2;
          atom_eelec[atom2] += delta2;
        // 2 - No reference, just cumulative Energy on atoms
        } else {
          // Cumulative evdw - divide between both atoms
          delta2 = e_vdw * 0.5;
          atom_evdw[atom1] += delta2;
          atom_evdw[atom2] += delta2;
          // Cumulative eelec - divide between both atoms
          delta2 = e_elec * 0.5;
          atom_eelec[atom1] += delta2;
          atom_eelec[atom2] += delta2;
        }
        Ncomparison++;
      } // End if not excluded
    } // End Inner loop
  } // End Outer loop
  //mprintf("\tPAIRWISE: ELJ = %12.4lf  Eelec = %12.4lf\n",ELJ,Eelec);
  //mprintf("\tPAIRWISE: %i interactions.\n",Ncomparison);

  // Print atoms for which the cumulative energy satisfies the given
  // cutoffs. Also create MOL2 files containing those atoms.
  AtomMask CutMask; // TEST
  std::vector<double> CutCharges; // TEST
  // EVDW
  if (Eout.IsOpen()) {
    if (RefFrame!=NULL)
      Eout.IO->Printf("\tPAIRWISE: Cumulative dEvdw:");
    else
      Eout.IO->Printf("\tPAIRWISE: Cumulative Evdw:");
    Eout.IO->Printf(" Evdw < %.4lf, Evdw > %.4lf\n",cut_evdw1,cut_evdw);
  }
  for (int atom = 0; atom < Parm->natom; atom++) {
    if (atom_evdw[atom]>cut_evdw || atom_evdw[atom]<cut_evdw1) {
      if (Eout.IsOpen()) 
        Eout.IO->Printf("\t\t%6i@%s: %12.4lf\n",atom+1,Parm->names[atom],atom_evdw[atom]);
      CutMask.AddAtom(atom);
      CutCharges.push_back(atom_evdw[atom]);
    }
  }
  if (cutout!=NULL && !CutMask.None()) {
    sprintf(buffer,"%s.evdw.mol2",cutout);
    if (WriteCutFrame(Parm, &CutMask, &CutCharges[0], frame, buffer)) return 1;
  }
  CutMask.Reset();
  CutCharges.clear();
  // EELEC
  if (Eout.IsOpen()) {
    if (RefFrame!=NULL)
      Eout.IO->Printf("\tPAIRWISE: Cumulative dEelec:");
    else
      Eout.IO->Printf("\tPAIRWISE: Cumulative Eelec:");
    Eout.IO->Printf(" Eelec < %.4lf, Eelec > %.4lf\n",cut_eelec1,cut_eelec);
  }
  for (int atom = 0; atom < Parm->natom; atom++) { 
    if (atom_eelec[atom]>cut_eelec || atom_eelec[atom]<cut_eelec1) {
      if (Eout.IsOpen()) 
        Eout.IO->Printf("\t\t%6i@%s: %12.4lf\n",atom+1,Parm->names[atom],atom_eelec[atom]);
      CutMask.AddAtom(atom);
      CutCharges.push_back(atom_eelec[atom]);
    }  
  }
  if (cutout!=NULL && !CutMask.None()) {
    sprintf(buffer,"%s.eelec.mol2",cutout);
    if (WriteCutFrame(Parm, &CutMask, &CutCharges[0], frame,buffer)) return 1;
  }
  // Reset cumulative energy arrays
  atom_eelec.assign(Parm->natom, 0);
  atom_evdw.assign(Parm->natom, 0);
  return 0;
}

/* Pairwise::action()
 */
int Pairwise::action() {
  if (Energy(&Mask0, F, P)) return 1;
  ds_vdw->Add(currentFrame, &ELJ);
  ds_elec->Add(currentFrame, &Eelec);

  return 0;
} 

/* Pairwise::print()
 */
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
