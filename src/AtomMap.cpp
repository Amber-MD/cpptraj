// AtomMap
#include <cstdlib>
#include <cstring>
#include "AtomMap.h"
#include "CpptrajStdio.h"

atommap::atommap() {
  M=NULL;
  natom=0;
  names=NULL;
  F=NULL;
  debug=1;
}

atommap::~atommap() {
  int i;
  if (M!=NULL) free(M);
  if (names!=NULL) {
    for (i=0; i<natom; i++) free(names[i]);
    free(names);
  }
}

/*------------------------------------------------
 * compareName(A1,A2,B1,B2) 
 * Compares pairs of names, return 1 if they match.  e.g. (A,B)==(A,B)==(B,A)
 */
int compareName(char *nameA1, char *nameA2, 
                const char *nameB1, const char *nameB2) {

  if ( ((strcmp(nameA1,nameB1)==0) && (strcmp(nameA2,nameB2)==0)) ||
       ((strcmp(nameA1,nameB2)==0) && (strcmp(nameA2,nameB1)==0)) )
    return 0;
  else
    return 1;
}
/*------------------------------------------------
 * getCut(atom1,atom2) 
 * Return a cutoff based on optimal covalent bond distance based on the 
 * identity of atom1 and atom2.
 * Treat X as chlorine for now.
 */
double atommap::getCut(char *atom1, char *atom2) {
  double cut;

  cut=1.60;

  /*self*/
  if (strcmp(atom1,atom2)==0) {
    if (strcmp(atom1,"H")==0) cut=0.74;
    if (strcmp(atom1,"N")==0) cut=1.45;
    if (strcmp(atom1,"C")==0) cut=1.54;
    if (strcmp(atom1,"O")==0) cut=0.74;
  }
  /* Others */
  else if ( compareName(atom1,atom2,"H","C")==0 )
    cut=1.09;
  else if ( compareName(atom1,atom2,"H","N")==0 )
    cut=1.01;
  else if ( compareName(atom1,atom2,"H","O")==0 )
    cut=0.96;
  else if ( compareName(atom1,atom2,"C","N")==0 )
    cut=1.47;
  else if ( compareName(atom1,atom2,"C","O")==0 )
    cut=1.43;
  else if ( compareName(atom1,atom2,"C","X")==0 )
    cut=1.76;
  else if ( compareName(atom1,atom2,"C","S")==0 )
    cut=1.83;
  else if ( compareName(atom1,atom2,"N","O")==0 )
    cut=1.47;
  else if ( compareName(atom1,atom2,"S","O")==0 )
    cut=1.48;
  else {
    if (debug>0) {
      mprintf("Warning: atommap::getCut: Cut not found for %s - %s\n",atom1,atom2);
      mprintf("                          Using default cutoff of %lf\n",cut);
    }
  }

  cut+=0.1;
  return cut;
}


/*------------------------------------------------
 * calcDist(M,natom)
 * Determine which atoms are bonded to each other in a given set of atoms
 * based on how close they are and their identity.
 */
int atommap::calcDist() {
  int i,j,nbondi,nbondj;
  double r,cut;

  for (i=0; i<natom-1; i++) {
    for(j=i+1; j<natom; j++) {
      if (debug>0) mprintf("%s_%i - %s_%i ",names[i],i,names[j],j);
      nbondi=M[i].nbond;
      nbondj=M[j].nbond;
      r=F->DIST(i,j);
      if (debug>0) mprintf("%lf ",r);
      /* Lookup bond distance based on atom names */
      cut=getCut(names[i],names[j]);
      if (r<cut) {
        if (debug>0) mprintf("nbondi=%i nbondj=%i ",nbondi,nbondj);
        if ((nbondi<4)&&(nbondj<4)) {
          M[i].bond[nbondi++]=j;
          M[j].bond[nbondj++]=i;
          M[i].nbond=nbondi;
          M[j].nbond=nbondj;
          if (debug>0) mprintf("BONDED!\n");
        } else
          if (debug>0) mprintf("MAXED!\n");
      } else
        if (debug>0) mprintf("NO_BOND!\n");
    } /* loop over j */
    if (debug>0) mprintf("\n"); 
  } /* loop over i */

  // DEBUG - print bonding information
  if (debug>0) {
    mprintf("atommap: Atom Bond information.\n");
    for (i=0; i<natom; i++) {
      mprintf("  Atom %s_%i has %i bonds.\n",names[i],i,M[i].nbond);
      for (j=0; j<M[i].nbond; j++) {
        nbondj=M[i].bond[j];
        mprintf("    to %s_%i\n",names[nbondj],nbondj);
      }
    }
  }

  return 0;
}

/*------------------------------------------------
 * compare(a,b)
 * Compare characters a and b, for use with qsort
 */
int compareChar (const void *a, const void *b) {

  return ( *(char*)a - *(char*)b );
}

/*------------------------------------------------
 * determineAtomID(M,natom)
 * Give each atom an identifier based on what atoms are bonded to it. The
 * first part is the atom itself, followed by an alphabetized list of 
 * bonded atoms. So C in O=C-H2 would be CHHO.
 */
void atommap::determineAtomID() {
  int i,j,atom;
  char *formula;

  formula=(char*) malloc(ATOMIDLENGTH*sizeof(char));
  if (debug>0) mprintf("ATOM IDs:\n");
  for (i=0; i<natom; i++) {
    strcpy(formula,"");
    for (j=0; j<M[i].nbond; j++) {
      atom=M[i].bond[j];
      strcat(formula,names[atom]);
    }
    qsort(formula,strlen(formula),sizeof(char),compareChar);
    strcpy(M[i].atomID,names[i]);
    strcat(M[i].atomID,formula);
    if (debug>0) mprintf("  Atom %i : %s\n",i,M[i].atomID);
  }
  free(formula);

  // Create a unique ID for each atom based on Atom IDs
  for (i=0; i<natom; i++) {
    memset(M[i].unique,' ',UNIQUELENGTH);
    strcpy(M[i].unique,M[i].atomID);
    for (j=0; j<M[i].nbond; j++) {
      atom=M[i].bond[j];
      strcat(M[i].unique,M[atom].atomID);
    }
    qsort(M[i].unique,strlen(M[i].unique),sizeof(char),compareChar);
  }

  // Determine which unique IDs are duplicated - set isUnique flag
  for (i=0; i<natom-1; i++) {
    for (j=i+1; j<natom; j++) {
      if (strcmp(M[i].unique,M[j].unique)==0) {
        // This unique string is duplicated, set isUnique flags
        M[i].isUnique=0;
        M[j].isUnique=0;
      }
    }
  }

  // Debug Output
  if (debug>0) {
    mprintf("UNIQUE IDs:\n");
    for (i=0; i<natom; i++) {
      mprintf("  Atom %i : %s",i,M[i].unique);
      if (M[i].isUnique) mprintf(" UNIQUE!");
      mprintf("\n");
    }
  }
  
}

/* atommap::setup()
 * Allocate memory for atom map. In order to easily create the uniqueID 
 * strings the atom names need to be 1 char long. Convert any multiple
 * char atom names to X for now.
 */
int atommap::setup() {
  int i;
  char *ptr;

  natom=P->natom;
  names=(char**) malloc(natom * sizeof(char*));
  // Set up atom names.
  for (i=0; i<natom; i++) {
    names[i]=(char*) malloc(2*sizeof(char));
    // position ptr at first non-space character in name
    ptr=P->names[i];
    while (*ptr==' ' && *ptr!='\0') ptr++;
    // if NULL something went wrong, abort
    if (*ptr=='\0') {
      strcpy(names[i],"");
      continue;
    }
    names[i][0]=ptr[0];
    // If C, check for L or l for chlorine
    if (ptr[0]=='C') {
      if (ptr[1]=='L' || ptr[1]=='l') names[i][0]='X';
    }
    names[i][1]='\0';
    // DEBUG
    mprintf("  Atom %i element: [%s]\n",i,names[i]);
  }
  // Allocate memory for atoms and initialize each atom
  // NOTE: visited still needed?
  M=(mapatom*) malloc( natom * sizeof(mapatom));
  for (i=0; i<natom; i++) {
    M[i].bond[0]=-1;
    M[i].bond[1]=-1;
    M[i].bond[2]=-1;
    M[i].bond[3]=-1;
    M[i].nbond=0;
    M[i].visited=0;
    memset(M[i].atomID,' ',ATOMIDLENGTH);
    memset(M[i].unique,' ',UNIQUELENGTH);
    M[i].isUnique=1; // Assume unique until proven otherwise
  }
  return 0;
}

// ============================================================================
/*
 * AtomMap::init()
 * Expected call: atommap <target> <reference> [out <filename>]
 */
int AtomMap::init() {
  char *refName, *targetName;
  int refIndex, targetIndex;
  int refatom,targetatom, *AtomMap,bond,numNonUniqueBond,nonUniqueID,j;

  targetName=A->getNextString();
  refName=A->getNextString();

  if (targetName==NULL) {
    mprintf("AtomMap::init: Error: No target specified.\n");
    return 1;
  }
  if (refName==NULL) {
    mprintf("AtomMap::init: Error: No reference specified.\n");
    return 1;
  }

  // Get reference and target frames
  refIndex=FL->GetFrameIndex(refName);
  //RefMap.F=FL->GetFrame(refName);
  RefMap.F=FL->GetFrame(refIndex);
  RefMap.P=FL->GetFrameParm(refIndex);
  if (RefMap.F==NULL || RefMap.P==NULL) {
    mprintf("AtomMap::init: Error: Could not get reference frame %s\n",refName);
    return 1;
  }
  targetIndex=FL->GetFrameIndex(targetName);
  //TargetMap.F=FL->GetFrame(targetName);
  TargetMap.F=FL->GetFrame(targetIndex);
  TargetMap.P=FL->GetFrameParm(targetIndex);
  if (TargetMap.F==NULL || TargetMap.P==NULL) {
    mprintf("AtomMap::init: Error: Could not get target frame %s\n",targetName);
    return 1;
  }

  RefMap.setup();
  RefMap.calcDist();
  RefMap.determineAtomID();

  TargetMap.setup();
  TargetMap.calcDist();
  TargetMap.determineAtomID();

  if (RefMap.natom!=TargetMap.natom) {
    mprintf("AtomMap::init: Warning: # atoms in reference %i not equal\n",RefMap.natom);
    mprintf("               to # atoms in target %i.\n",TargetMap.natom);
  }

  /* Atoms have now been assigned IDs. Match up the unique strings in Ref with 
   * unique strings in target.
   */
  AtomMap=(int*) malloc( RefMap.natom * sizeof(int));
  for (refatom=0; refatom<RefMap.natom; refatom++) {
    AtomMap[refatom]=-1;
    if (RefMap.M[refatom].isUnique) {
      for (targetatom=0; targetatom<TargetMap.natom; targetatom++) {
        if (TargetMap.M[targetatom].isUnique) {
          if ( strcmp(TargetMap.M[targetatom].unique,
                      RefMap.M[refatom].unique)==0 ) {
            AtomMap[refatom]=targetatom;
            mprintf("* TargetAtom %i(%s) maps to RefAtom %i(%s)\n",
                    targetatom,TargetMap.P->names[targetatom],
                    refatom,RefMap.P->names[refatom]);
          } // If unique strings match
        } // If target atom is unique
      } // Loop over target atoms
    } // If reference atom is unique
  } // Loop over reference atoms

  /* Two potential strategies for assigning the remainder: (1) try to RMS
   * fit the unique atoms, then assign based on distances; or (2) try to 
   * assign each non-unique atom based on bonding information to unique
   * atoms.
   */

  /* Search for unique atoms that have only one non-unique bound.
   */
  for (refatom=0; refatom<RefMap.natom; refatom++) {
    numNonUniqueBond=0; 
    nonUniqueID=-1;
    if (RefMap.M[refatom].isUnique) {
      for (bond=0; bond<RefMap.M[refatom].nbond; bond++) {
        targetatom=RefMap.M[refatom].bond[bond];
        if (RefMap.M[targetatom].isUnique==0) {
          nonUniqueID=targetatom;
          numNonUniqueBond++;
        }
      }
      if (numNonUniqueBond==1) {
        mprintf("  Nonunique RefAtom %i(%s) bound to unique RefAtom %i(%s).\n",
                nonUniqueID,RefMap.P->names[nonUniqueID],
                refatom,RefMap.P->names[refatom]);
        // Find the same unique atom in Target
        targetatom=AtomMap[refatom];
        /* Find the target atom bonded to unique target atom with same unique ID
         * as reference atom bonded to unique reference atom.
         */
        for (bond=0; bond<TargetMap.M[targetatom].nbond; bond++) {
           j=TargetMap.M[targetatom].bond[bond];
           if (strcmp(RefMap.M[nonUniqueID].unique,
                      TargetMap.M[j].unique)==0     ) {
             AtomMap[nonUniqueID]=j;
           }
        }
      }
    }
  }
  
  // DEBUG
  for (refatom=0; refatom<RefMap.natom; refatom++) {
    targetatom=AtomMap[refatom];
    if (targetatom>=0)
      mprintf("* TargetAtom %i(%s) maps to RefAtom %i(%s)\n",
                      targetatom,TargetMap.P->names[targetatom],
                      refatom,RefMap.P->names[refatom]);
  }        

  free(AtomMap);

  return 0;
}

