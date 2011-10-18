// Pairwise
#include <cmath>
#include "Action_Pairwise.h"
#include "CpptrajStdio.h"
#include "vectormath.h"

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
} 

// DESTRUCTOR
Pairwise::~Pairwise() {
  //fprintf(stderr,"Pairwise Destructor.\n");
  if (skipv!=NULL) delete[] skipv;
  if (natexidx!=NULL) delete[] natexidx;
}

/* Pairwise::init()
 * Expected call: pairwise [<mask>] eout <filename>
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int Pairwise::init( ) {
  char *mask0, *eout, *refmask, *referenceName;
  int refindex;

  // Get Keywords
  eout = A->getKeyString("eout",NULL);
  if (eout==NULL) {
    mprinterr("Error: Pairwise::init(): No output file specified (use 'eout' keyword).\n");
    return 1;
  }
  referenceName=A->getKeyString("ref",NULL);
  refindex=A->getKeyInt("refindex",-1);
  
  // Get Masks
  mask0 = A->getNextMask();
  //fprintf(stdout,"    Mask 1: %s\n",mask0);
  Mask0.SetMaskString(mask0);
  refmask = A->getNextMask();
  if (refmask!=NULL)
    RefMask.SetMaskString(refmask);
  else
    RefMask.SetMaskString(mask0);

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
    if (RefMask.SetupMask(RefParm, debug)) return 1;
    if (RefMask.None()) {
      mprinterr("    Error: Pairwise::init: No atoms selected in reference mask.\n");
      return 1;
    }
    // Allocate exclusion list for reference
    AllocateExclusion(RefParm);
    // Calculate energy for reference
    Energy(&RefMask, RefFrame, RefParm);
    mprintf("DEBUG:\tReference ELJ= %12.4lf  Eelec= %12.4lf\n",ELJ,Eelec); 
  }

  // Datasets
  ds_vdw = Eout.Add(DOUBLE, (char*)NULL, "EVDW");
  ds_elec= Eout.Add(DOUBLE, (char*)NULL,"EELEC");
  // Add datasets to data file list
  DFL->Add(eout,ds_vdw);
  DFL->Add(eout,ds_elec);

  mprintf("    PAIRWISE: Atoms in mask %s, output to %s.\n",Mask0.maskString,eout);

  return 0;
}

/* Pairwise::AllocateExclusion()
 * Allocate memory for the exclusion list based on the given parm.
 */
int Pairwise::AllocateExclusion(AmberParm *Parm) {
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
    int idx = 0;
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

/* Pairwise::setup()
 * Set up mask, allocate memory for exclusion list.
 */
int Pairwise::setup() {
  // Set up mask
  if ( Mask0.SetupMask(P,debug) ) return 1;
  if (Mask0.None()) {
    mprintf("    Error: Pairwise::setup: Mask has no atoms.\n");
    return 1;
  }
  // Allocate memory for exclusion list
  if (AllocateExclusion(P)) return 1; 
  // Print info for this parm
  mprintf("    PAIRWISE: Mask %s corresponds to %i atoms.\n",Mask0.maskString, Mask0.Nselected);
        
  return 0;  
}

/* Pairwise::SetupExclusion()
 * Sets up the exclusion list, skipv, for the given atom; skipv will be
 * true for atoms that will be skipped (i.e. excluded) in the non-bonded
 * calc. skipv and natexidx should be set up prior to calling this routine.
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
  Parm->GetLJparam(&Acoef, &Bcoef,atom1,atom2);
  r2=1/rij2;
  r6=r2*r2*r2;
  r12=r6*r6;
  f12=Acoef*r12;               // A/r^12
  f6=Bcoef*r6;                 // B/r^6
  energy=f12-f6;               // (A/r^12)-(B/r^6)
  //ELJ+=energy;
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
  // NOTE: Currently cpptraj converts charge to units of e- when topology
  // is read in. Need to convert back here. Not really efficient, just a 
  // hack for now.
  qi = Parm->charge[atom1] * 18.2223;
  qj = Parm->charge[atom2] * 18.2223;
  qiqj = qi * qj;
  energy=kes * (qiqj/rij);
  //Eelec+=energy;
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

/* Pairwise::Energy()
 * Calculate pairwise energy for the given mask, frame, and parm. Sets
 * ELJ and Eelec.
 */
void Pairwise::Energy(AtomMask *atommask, Frame *frame, AmberParm *Parm) {
  int atom1, atom2;
  double rij, rij2, JI[3], force;
  // Loop over all atom pairs excluding self
  ELJ = 0;
  Eelec = 0;
  JI[0]=0; JI[1]=0; JI[2]=0;
  // Outer loop
  for (int maskidx1 = 0; maskidx1 < atommask->Nselected - 1; maskidx1++) {
    atom1 = atommask->Selected[maskidx1];
    if (debug>0) mprintf("\tPAIRWISE: ATOM %i\n",atom1);
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
        ELJ += Energy_LJ(Parm, atom1, atom2, rij2, &force);
        // Coulomb Energy
        Eelec += Energy_Coulomb(Parm, atom1, atom2, rij, &force);
      }
    } // End Inner loop

  } // End Outer loop
  //mprintf("\tPAIRWISE: ELJ = %12.4lf  Eelec = %12.4lf\n",ELJ,Eelec);
}

/* Pairwise::action()
 */
int Pairwise::action() {
  Energy(&Mask0, F, P);
  ds_vdw->Add(currentFrame, &ELJ);
  ds_elec->Add(currentFrame, &Eelec);
  
  return 0;
} 

