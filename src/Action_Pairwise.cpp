// Pairwise
#include <cmath>
#include "Action_Pairwise.h"
#include "CpptrajStdio.h"
#include "vectormath.h"

// CONSTRUCTOR
Pairwise::Pairwise() {
  //fprintf(stderr,"Pairwise Con\n");
  skipv = NULL;
  natexidx = NULL;
  hasExclusion=true;
  kes = 1.0;
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
  char *mask0, *eout;

  // Get Keywords
  eout = A->getKeyString("eout",NULL);
  if (eout==NULL) {
    mprinterr("Error: Pairwise::init(): No output file specified (use 'eout' keyword).\n");
    return 1;
  }
  
  // Get Masks
  mask0 = A->getNextMask();
  //fprintf(stdout,"    Mask 1: %s\n",mask0);
  Mask0.SetMaskString(mask0);

  // Datasets
  ds_vdw = Eout.Add(DOUBLE, (char*)NULL, "EVDW");
  ds_elec= Eout.Add(DOUBLE, (char*)NULL,"EELEC");
  // Add datasets to data file list
  DFL->Add(eout,ds_vdw);
  DFL->Add(eout,ds_elec);

  mprintf("    PAIRWISE: Atoms in mask %s, output to %s.\n",Mask0.maskString,eout);

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
  if (skipv!=NULL) delete[] skipv;
  skipv = new bool[ P->natom ];
  for (int atom=0; atom < P->natom; atom++) skipv[atom]=false;
  // Check if exclusion info present
  if (P->NumExcludedAtoms(0)==-1) hasExclusion=false;
  if (hasExclusion) {
    // Create an array holding indices for each atom into NATEX
    natexidx = new int[ P->natom ];
    int idx = 0;
    for (int atom=0; atom < P->natom; atom++) {
      natexidx[atom] = idx;
      idx += P->NumExcludedAtoms(atom);
    }
    // DEBUG
    if (debug>0) {
      mprintf("DEBUG: EXCLUDED ATOM IDX:\n");
      for (int atom=0; atom < P->natom; atom++) 
        mprintf("\tAtom %i: Index %i\n",atom,natexidx[atom]);
    }
  }
  // Check if LJ parameters present - need at least 2 atoms for it to matter.
  if (P->natom>1) {
    double Atemp = 0;
    double Btemp = 0;
    if (P->GetLJparam(&Atemp, &Btemp, 0, 1)) {
      mprinterr("Error: Pairwise::setup(): Parm does not have LJ information.\n");
      return 1;
    }
  }
  
  // Print info for this parm
  mprintf("    PAIRWISE: Mask %s corresponds to %i atoms.\n",Mask0.maskString, Mask0.Nselected);
        
  return 0;  
}

/* Pairwise::action()
 */
int Pairwise::action() {
  int atom1,atom2;
  double rij, rij2, JI[3];
  double energy, force;
  double Acoef, Bcoef, ELJ;
  double Eelec;
  // Loop over all atom pairs excluding self
  ELJ = 0;
  Eelec = 0;
  JI[0]=0; JI[1]=0; JI[2]=0;
  // Outer loop
  for (int maskidx1 = 0; maskidx1 < Mask0.Nselected - 1; maskidx1++) {
    atom1 = Mask0.Selected[maskidx1];
    if (debug>0) mprintf("\tPAIRWISE: ATOM %i\n",atom1);
    if (hasExclusion) {
      // Set up exclusion list for atom1
      for (int atom2=atom1+1; atom2 < P->natom; atom2++) skipv[atom2]=false;
      int jexcl = natexidx[atom1];
      int jexcl_last = jexcl + P->NumExcludedAtoms(atom1);
      for (int idx = jexcl; idx < jexcl_last; idx++) {
        int natex = P->Natex(idx);
        if (debug>0) mprintf("\t\tPair %i - %i will be excluded.\n",atom1,natex);
        //if (natex < 0) {
        //  mprinterr("Error getting excluded atom index %i for atom %i\n",idx,atom1);
        //  return 1;
        //}
        if (natex>=0) skipv[ natex ] = true;
      }
    }
    // Inner loop
    for (int maskidx2 = maskidx1 + 1; maskidx2 < Mask0.Nselected; maskidx2++) {
      atom2 = Mask0.Selected[maskidx2];
      if (debug>0) mprintf("\t\tCalc for Pair %i - %i:\n",atom1,atom2);
      int coord1 = atom1 * 3;
      int coord2 = atom2 * 3;
      // Calculate the vector pointing from atom2 to atom1
      vector_sub(JI, F->X+coord1, F->X+coord2);
      // Normalize
      rij = vector_norm(JI, &rij2);
      // Non-bonded energy
      if (!skipv[atom2]) {
        // Lennard Jones 6-12 Energy
        P->GetLJparam(&Acoef, &Bcoef,atom1,atom2);
        double r2=1/rij2;
        double r6=r2*r2*r2;
        double r12=r6*r6;
        double f12=Acoef*r12; // A/r^12
        double f6=Bcoef*r6;   // B/r^6
        energy=f12-f6;        // (A/r^12)-(B/r^6)
        ELJ+=energy;
        // LJ Force 
        force=((12*f12)-(6*f6))*r2; // (12A/r^13)-(6B/r^7)
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
        // Coulomb Energy
        // NOTE: Currently cpptraj converts charge to units of e- when topology
        // is read in. Need to convert back here. Not really efficient, just a 
        // hack for now.
        double qi = P->charge[atom1] * 18.2223;
        double qj = P->charge[atom2] * 18.2223;
        double qiqj = qi * qj;
        energy=kes * (qiqj/rij);
        Eelec+=energy;
        // Coulomb Force
        force=energy/rij; // kes*(qiqj/r)*(1/r)
        //scalarmult(f,JI,F);
        //M[i].f[0]+=f[0];
        //M[i].f[1]+=f[1];
        //M[i].f[2]+=f[2];
        //M[j].f[0]-=f[0];
        //M[j].f[1]-=f[1];
        //M[j].f[2]-=f[2];
        if (debug>0) {
          mprintf("\t\t\tCoulomb Force:  q%i= %lf  q%i= %lf  F= %lf  E= %lf\n",
                  atom1,P->charge[atom1],atom2,P->charge[atom2],force,energy);
          //mprintf("\t\t\tAtom %i : Fx= %lf  Fy= %lf  Fz= %lf\n",atom1,f[0],f[1],f[2]);
          //mprintf("\t\t\tAtom %i : Fx= %lf  Fy= %lf  Fz= %lf\n",atom2,-f[0],-f[1],-f[2]);
        }
      }
    } // End Inner loop

  } // End Outer loop
  //mprintf("\tPAIRWISE: ELJ = %12.4lf  Eelec = %12.4lf\n",ELJ,Eelec);
  ds_vdw->Add(currentFrame, &ELJ);
  ds_elec->Add(currentFrame, &Eelec);
  
  return 0;
} 

