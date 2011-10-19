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
  isReference=false;
} 

// DESTRUCTOR
Pairwise::~Pairwise() {
  //fprintf(stderr,"Pairwise Destructor.\n");
  if (skipv!=NULL) delete[] skipv;
  if (natexidx!=NULL) delete[] natexidx;
}

/* Pairwise::init()
 * Expected call: pairwise [<name>] [<mask>] [out <filename>] [cuteelec <cute>] [cutevdw <cutv>]
                           [ref <reffilename> | refindex <ref#>]
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int Pairwise::init( ) {
  char *mask0, *eout, *refmask, *referenceName, *ds_name;
  int refindex;

  // Get Keywords
  eout = A->getKeyString("out",NULL);
  referenceName=A->getKeyString("ref",NULL);
  refindex=A->getKeyInt("refindex",-1);
  cut_eelec = A->getKeyDouble("cuteelec",1.0);
  cut_eelec1 = -cut_eelec;
  cut_evdw = A->getKeyDouble("cutevdw",1.0);
  cut_evdw1 = -cut_evdw;
  
  // Get Masks
  mask0 = A->getNextMask();
  //fprintf(stdout,"    Mask 1: %s\n",mask0);
  Mask0.SetMaskString(mask0);
  refmask = A->getNextMask();
  if (refmask!=NULL)
    RefMask.SetMaskString(refmask);
  else
    RefMask.SetMaskString(mask0);

   // Datasets
  ds_name = A->getNextString();
  ds_vdw = DSL->AddMulti(DOUBLE, ds_name, "EVDW");
  ds_elec= DSL->AddMulti(DOUBLE, ds_name, "EELEC");
  // Add datasets to data file list
  DFL->Add(eout,ds_vdw);
  DFL->Add(eout,ds_elec);

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
    // Set up comparison array
    int N_interactions = NumInteractions(&RefMask, RefParm);
    ref_eelec.clear();
    ref_eelec.resize( N_interactions, 0 );
    ref_evdw.clear();
    ref_evdw.resize( N_interactions, 0 );
    // Calculate energy for reference
    // Setting isReference to true lets energy routines know to fill the
    // ref_X arrays.
    isReference = true;
    Energy(&RefMask, RefFrame, RefParm);
    mprintf("DEBUG:\tReference ELJ= %12.4lf  Eelec= %12.4lf\n",ELJ,Eelec);
    mprintf("\tSize of reference eelec array: %u\n",ref_eelec.size());
    mprintf("\tSize of reference evdw array: %u\n",ref_evdw.size());
    isReference=false;
  }

  mprintf("    PAIRWISE: Atoms in mask %s, output to %s.\n",Mask0.maskString,eout);
  if (RefFrame!=NULL) {
    mprintf("              Reference index %i\n",refindex);
    mprintf("              Eelec absolute cutoff: %12.4lf\n",cut_eelec);
    mprintf("              Evdw absolute cutoff: %12.4lf\n",cut_evdw);
  }

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
  if ( Mask0.SetupMask(P,debug) ) return 1;
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
    if (N_interactions != (int) ref_eelec.size()) {
      mprinterr("Error: Pairwise::action: Size of ref %u != size of frame %i.\n",
                ref_eelec.size(), N_interactions);
      return 1;
    }
    atom_eelec.clear();
    atom_eelec.resize(P->natom, 0);
    atom_evdw.clear();
    atom_evdw.resize(P->natom, 0);
  }
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
  Parm->GetLJparam(&Acoef, &Bcoef,atom1,atom2);
  r2=1/rij2;
  r6=r2*r2*r2;
  r12=r6*r6;
  f12=Acoef*r12;               // A/r^12
  f6=Bcoef*r6;                 // B/r^6
  energy=f12-f6;               // (A/r^12)-(B/r^6)
  // Calculate difference to reference
  if (RefFrame!=NULL) {
    if (!isReference) {
      double delta = ref_evdw[ Ncomparison ] - energy;
      if (delta > cut_evdw || delta < cut_evdw1) {
        mprintf("\tAtom %6i@%4s-%6i@%4s dEvdw= %12.4lf\n",atom1+1,Parm->names[atom1],
                atom2+1,Parm->names[atom2],delta);
      }
      // Divide the total pair energy between both atoms.
      double delta2 = delta * 0.5;
      atom_evdw[atom1] += delta2;
      atom_evdw[atom2] += delta2;
    } else {
      ref_evdw[ Ncomparison ] = energy;
    }
  }
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
  // Calcualte difference to reference
  if (RefFrame!=NULL) {
    if (!isReference) {
      double delta = ref_eelec[ Ncomparison ] - energy;
      if (delta > cut_eelec || delta < cut_eelec1) {
        mprintf("\tAtom %6i@%4s-%6i@%4s dEelec= %12.4lf\n",atom1+1,Parm->names[atom1],
                atom2+1,Parm->names[atom2],delta);
      }
      // Divide the total pair energy between both atoms.
      double delta2 = delta * 0.5;
      atom_eelec[atom1] += delta2;
      atom_eelec[atom2] += delta2;
    } else {
      ref_eelec[ Ncomparison ] = energy;
    }
  }
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

/* Pairwise::WriteCutFrame()
 */
void Pairwise::WriteCutFrame(AmberParm *Parm, AtomMask *CutMask, double *CutCharges,
                             Frame *frame, char *outfilename) 
{
    AmberParm *CutParm;
    Frame *CutFrame;
    TrajectoryFile tout;
    // TEST: Write file containing only cut atoms
    CutParm = Parm->modifyStateByMask(CutMask->Selected, CutMask->Nselected);
    CutParm->SetCharges(CutCharges);
    CutFrame = new Frame(CutParm->natom, CutParm->mass);
    CutFrame->SetFrameFromMask(frame, CutMask);
    tout.SetupWrite(outfilename,NULL,CutParm,MOL2FILE);
    tout.WriteFrame(0,CutParm,CutFrame->X, NULL, NULL, 0.0);
    tout.EndTraj();
    delete CutParm;
    delete CutFrame;
}

/* Pairwise::Energy()
 * Calculate pairwise energy for the given mask, frame, and parm. Sets
 * ELJ and Eelec.
 */
void Pairwise::Energy(AtomMask *atommask, Frame *frame, AmberParm *Parm) {
  int atom1, atom2;
  double rij, rij2, JI[3], force;
  char buffer[256]; // NOTE: Temporary
  // Loop over all atom pairs excluding self
  ELJ = 0;
  Eelec = 0;
  JI[0]=0; JI[1]=0; JI[2]=0;
  Ncomparison = 0;
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
        // DEBUG - safety valve
        if (RefFrame != NULL && Ncomparison >= (int) ref_eelec.size()) {
          mprinterr("Internal Error: Attempting to access element %i of array with size %u\n",
                    Ncomparison, ref_eelec.size());
          mprintf("Internal Error: Attempting to access element %i of array with size %u\n",
                  Ncomparison, ref_eelec.size());
          mprintf("atom1: %i  atom2: %i\n",atom1,atom2);
          return;
        }
        Eelec += Energy_Coulomb(Parm, atom1, atom2, rij, &force);
        Ncomparison++;
      }
    } // End Inner loop

  } // End Outer loop
  //mprintf("\tPAIRWISE: ELJ = %12.4lf  Eelec = %12.4lf\n",ELJ,Eelec);
  //mprintf("\tPAIRWISE: %i interactions.\n",Ncomparison);

  if (!isReference && RefFrame!=NULL) {
    AtomMask CutMask; // TEST
    std::vector<double> CutCharges; // TEST

    mprintf("\tPAIRWISE: Cumulative dEvdw:\n");
    for (int atom = 0; atom < Parm->natom; atom++) {
      if (atom_evdw[atom]>cut_evdw || atom_evdw[atom]<cut_evdw1) {
        mprintf("\t\t%6i@%s: %12.4lf\n",atom,Parm->names[atom],atom_evdw[atom]);
        CutMask.AddAtom(atom);
        CutCharges.push_back(atom_evdw[atom]);
      }
    }
    sprintf(buffer,"evdw.test.mol2");
    WriteCutFrame(Parm, &CutMask, &CutCharges[0], frame, buffer);
    CutMask.Reset();
    CutCharges.clear();

    mprintf("\tPAIRWISE: Cumulative dEelec:\n");
    for (int atom = 0; atom < Parm->natom; atom++) { 
      if (atom_eelec[atom]>cut_eelec || atom_eelec[atom]<cut_eelec1) {
        mprintf("\t\t%6i@%s: %12.4lf\n",atom,Parm->names[atom],atom_eelec[atom]);
        // TEST: Create mask for keeping only atoms that satisfy cut
        CutMask.AddAtom(atom);
        CutCharges.push_back(atom_eelec[atom]);
      }  
    }
    // TEST: Write file containing only cut atoms
    sprintf(buffer,"eelec.test.mol2");
    WriteCutFrame(Parm, &CutMask, &CutCharges[0], frame,buffer);
/*    AmberParm *CutParm = P->modifyStateByMask(CutMask.Selected, CutMask.Nselected);
    CutParm->SetCharges(&CutCharges[0]);
    Frame *CutFrame = new Frame(CutParm->natom, CutParm->mass);
    CutFrame->SetFrameFromMask(frame, &CutMask);
    TrajectoryFile tout;
    tout.SetupWrite((char*)"test.mol2",NULL,CutParm,MOL2FILE);
    tout.WriteFrame(0,CutParm,CutFrame->X, NULL, NULL, 0.0);
    tout.EndTraj();
    delete CutParm;
    delete CutFrame;*/
    // Reset
    atom_eelec.assign(Parm->natom, 0);
    atom_evdw.assign(Parm->natom, 0);
  }

}

/* Pairwise::action()
 */
int Pairwise::action() {
  Energy(&Mask0, F, P);
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
