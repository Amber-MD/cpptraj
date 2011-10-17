// Pairwise
#include <cmath>
#include "Action_Pairwise.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Pairwise::Pairwise() {
  //fprintf(stderr,"Pairwise Con\n");
  skipv = NULL;
  natexidx = NULL;
  hasExclusion=true;
} 

// DESTRUCTOR
Pairwise::~Pairwise() {
  //fprintf(stderr,"Pairwise Destructor.\n");
  if (skipv!=NULL) delete[] skipv;
  if (natexidx!=NULL) delete[] natexidx;
}

/* Pairwise::init()
 * Expected call: pairwise <mask> <mask2> [out filename] [geom] [noimage]
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int Pairwise::init( ) {
  char *mask0;

  // Get Keywords
  
  // Get Masks
  mask0 = A->getNextMask();
  //fprintf(stdout,"    Mask 1: %s\n",mask0);
  Mask0.SetMaskString(mask0);

  // Dataset 
  //dist = DSL->Add(DOUBLE, A->getNextString(),"Dis");
  //if (dist==NULL) return 1;
  // Add dataset to data file list
  //DFL->Add(distanceFile,dist);

  mprintf("    PAIRWISE: Atoms in mask %s\n",Mask0.maskString);

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
    mprintf("DEBUG: EXCLUDED ATOM IDX:\n");
    for (int atom=0; atom < P->natom; atom++) 
      mprintf("\tAtom %i: Index %i\n",atom,natexidx[atom]);
  }

  // Print info for this parm
  mprintf("    PAIRWISE: Mask %s corresponds to %i atoms.\n",Mask0.maskString, Mask0.Nselected);
        
  return 0;  
}

/* Pairwise::action()
 */
int Pairwise::action() {
  int atom1;
  // Loop over all atom pairs excluding self
  for (int maskidx1 = 0; maskidx1 < Mask0.Nselected - 1; maskidx1++) {
    atom1 = Mask0.Selected[maskidx1];
    mprintf("\tPAIRWISE: ATOM %i\n",atom1);
    if (hasExclusion) {
      // Set up exclusion list for atom1
      for (int atom2=atom1+1; atom2 < P->natom; atom2++) skipv[atom2]=false;
      int jexcl = natexidx[atom1];
      int jexcl_last = jexcl + P->NumExcludedAtoms(atom1);
      for (int idx = jexcl; idx < jexcl_last; idx++) {
        int natex = P->Natex(idx);
        mprintf("\t\tPair %i - %i will be excluded.\n",atom1,natex);
        if (natex < 0) {
          mprinterr("Error getting excluded atom index %i for atom %i\n",idx,atom1);
          return 1;
        }
        skipv[ natex ] = true;
      }
    }

  }
  
  return 0;
} 

