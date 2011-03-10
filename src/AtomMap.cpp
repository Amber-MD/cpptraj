// AtomMap
#include <cstdlib>
#include <cstring>
#include "AtomMap.h"
#include "CpptrajStdio.h"
// DEBUG - pdb write
#include "PDBfileRoutines.h"

// atommap CONSTRUCTOR
atommap::atommap() {
  M=NULL;
  natom=0;
  names=NULL;
  F=NULL;
  debug=1;
}

// atommap DESTRUCTOR
atommap::~atommap() {
  int i;
  if (M!=NULL) free(M);
  if (names!=NULL) {
    for (i=0; i<natom; i++) free(names[i]);
    free(names);
  }
}

//--------- PRIVATE ROUTINES ---------------------------------------
/*
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

/*
 * compare(a,b)
 * Compare characters a and b, for use with qsort
 */
int compareChar(const void *a, const void *b) {
  return ( *(char*)a - *(char*)b );
}
//------------------------------------------------------------------

/*
 * atommap::getCut() 
 * Return a cutoff based on optimal covalent bond distance based on the 
 * identities of atom1 and atom2.
 * Treat X as chlorine for now.
 */
double atommap::getCut(char *atom1, char *atom2) {
  double cut;

  cut=1.60;

  // Self
  if (strcmp(atom1,atom2)==0) {
    if (strcmp(atom1,"H")==0) cut=0.74;
    if (strcmp(atom1,"N")==0) cut=1.45;
    if (strcmp(atom1,"C")==0) cut=1.54;
    if (strcmp(atom1,"O")==0) cut=0.74;
  }
  // Others 
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

/*
 * atommap::calcDist()
 * Determine which atoms are bonded to each other in a given set of atoms
 * based on how close they are and their identity.
 */
int atommap::calcDist() {
  int i,j,nbondi,nbondj;
  double r,cut;

  for (i=0; i<natom-1; i++) {
    for(j=i+1; j<natom; j++) {
      if (debug>1) mprintf("%s_%i - %s_%i ",names[i],i,names[j],j);
      nbondi=M[i].nbond;
      nbondj=M[j].nbond;
      r=F->DIST(i,j);
      if (debug>1) mprintf("%lf ",r);
      // Lookup bond distance based on atom names 
      cut=getCut(names[i],names[j]);
      if (r<cut) {
        if (debug>1) mprintf("nbondi=%i nbondj=%i ",nbondi,nbondj);
        if ((nbondi<4)&&(nbondj<4)) {
          M[i].bond[nbondi++]=j;
          M[j].bond[nbondj++]=i;
          M[i].nbond=nbondi;
          M[j].nbond=nbondj;
          if (debug>1) mprintf("BONDED!\n");
        } else
          if (debug>1) mprintf("MAXED!\n");
      } else
        if (debug>1) mprintf("NO_BOND!\n");
    } // END loop over j
    if (debug>1) mprintf("\n"); 
  } // END loop over i

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

/*
 * atommap::printBonds()
 * Print bond information for each atom in the map, indicate whether
 * atoms are unique or not.
 */
void atommap::printBonds() {
  int atom,bond,bondatom;

  for (atom=0; atom<natom; atom++) {
    mprintf("  Atom %4i: %s-%i |",atom,names[atom],M[atom].isUnique);
    for (bond=0; bond<M[atom].nbond; bond++) {
      bondatom = M[atom].bond[bond];
      mprintf(" %s-%i",names[bondatom],M[bondatom].isUnique);
    }
    mprintf("\n");
  }
}

/*
 * atommap::determineAtomID()
 * Give each atom an identifier based on what atoms are bonded to it. The
 * first part is the atom itself, followed by an alphabetized list of 
 * bonded atoms. So C in O=C-H2 would be CHHO.
 * Then determine which identifier strings are unique. 
 */
void atommap::determineAtomID() {
  int i,j,atom;
  char *formula;
  bool isRepeated;
  int k, atom2;

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
        // This unique string is duplicated, set isUnique flags to false
        M[i].isUnique=0;
        M[j].isUnique=0;
      }
    }
  }

  // For each atom with a truly unique ID, determine if it is bonded to a
  // non-unique partner. If that partner is itself unique among bonded
  // partners (e.g. H2-C-N where C is unique, N is unique by extension),
  // give it a unique ID of atomID-element
  for (i=0; i<natom; i++) {
    if (M[i].isUnique) {
      // Check bonds of unique atom i for non-unique
      for (j=0; j<M[i].nbond; j++) {
        atom = M[i].bond[j];
        if (!M[atom].isUnique) {
          // Check if non-unique atom name is same as other atoms bonded to atom i
          isRepeated=false;
          for (k=0; k<M[i].nbond; k++) {
            atom2 = M[i].bond[k];
            if (atom==atom2) continue;
            if (M[atom2].isUnique) continue;
            if ( strcmp(names[atom],names[atom2])==0 ) {
              isRepeated=true;
              break;
            }
          } // END loop k over bonds of atom i
          // If non-unique atom is not repeated, give it a unique ID
          if (!isRepeated) {
            strcpy(M[atom].unique,M[i].unique);
            strcat(M[atom].unique,"-");
            strcat(M[atom].unique,names[atom]);
            M[atom].isUnique=1;
          }
        } // End if bonded atom is not unique
      } // End loop j over bonds of atom i
    } // End if atom i is unique
  } // End loop i over atoms in map 

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

/* 
 * atommap::setup()
 * Allocate memory for atom map. In order to easily create the uniqueID 
 * strings the atom names need to be 1 char long. Convert chlorine 
 * to X for now.
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

/*
 * atommap::fourAtoms()
 * Fill the passed in int array with indices of four unique atoms. The first 
 * atom should be bonded to at least 1 non-unique atom. If all atoms are
 * unique, return 1, otherwise return 0;
 */
int atommap::fourAtoms(int *Atom) {
  int n, atom, bond, bondedatom;

  if (Atom==NULL) return 2;
  Atom[0]=-1; Atom[1]=-1; Atom[2]=-1; Atom[3]=-1;

  n=0;
  for (n=0; n<4; n++) {
    Atom[n]=-1;
    // If this is first atom, find unique atom bonded to non-unique atom
    // If not, make sure this atom is unique and bonded to Atom[n-1]
    for (atom=0; atom<natom; atom++) {
      if (!M[atom].isUnique) continue;
      for (bond=0; bond<M[atom].nbond; bond++) {
        bondedatom=M[atom].bond[bond];

        if (n==0) {
          if (!M[bondedatom].isUnique) {
            Atom[0]=atom;
            break;
          }
        } else {
          if (bondedatom==Atom[n-1]) {
            Atom[n]=atom;
            break;
          }
        }

      }
      // If atom has been found leave loop over natom
      if (Atom[n]!=-1) break;
    }
  }
  mprintf("    First unique atom bonded to a non-unique atom: %i\n",Atom[0]);
  if (Atom[0]==-1) return 1;

  // Find 4 unique atoms bonded to first atom


  return 0;
}

// ============================================================================
/*
 * AtomMap::init()
 * Expected call: atommap <target> <reference>
 * Attempt to create a map from atoms in target to atoms in reference solely
 * based on how they are bonded (not how they are named). 
 */
int AtomMap::init() {
  char *refName, *targetName;
  int refIndex, targetIndex;
  int refatom,targetatom, *AtomMap;
  // For RMSD
  int fouratomlist[4];
  int numMappedAtoms=0;
  Frame *refFrame, *tgtFrame;
  double rms, Rot[9],Trans[6];
  // DEBUG - pdb write
  char buffer[83];

  // Get Args
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

  // Get reference index based on filename 
  refIndex=FL->GetFrameIndex(refName);
  // Get reference frame
  RefMap.F=FL->GetFrame(refIndex);
  // Get reference parm
  RefMap.P=FL->GetFrameParm(refIndex);
  if (RefMap.F==NULL || RefMap.P==NULL) {
    mprintf("AtomMap::init: Error: Could not get reference frame %s\n",refName);
    return 1;
  }
  // Get target index based on filename
  targetIndex=FL->GetFrameIndex(targetName);
  // Get target frame 
  TargetMap.F=FL->GetFrame(targetIndex);
  // Get target parm
  TargetMap.P=FL->GetFrameParm(targetIndex);
  if (TargetMap.F==NULL || TargetMap.P==NULL) {
    mprintf("AtomMap::init: Error: Could not get target frame %s\n",targetName);
    return 1;
  }

  // For each map, set up (get element for each atom, initialize map mem),
  // determine what atoms are bonded to each other via simple distance
  // cutoffs, the give each atom an ID based on what atoms are bonded to
  // it, noting which IDs are unique for that map.

  RefMap.setup();
  RefMap.calcDist();
  RefMap.determineAtomID();
  RefMap.printBonds();

  TargetMap.setup();
  TargetMap.calcDist();
  TargetMap.determineAtomID();
  TargetMap.printBonds();

  // Number of atoms in each map MUST be equal
  if (RefMap.natom!=TargetMap.natom) {
    mprintf("AtomMap::init: Warning: # atoms in reference %i not equal\n",RefMap.natom);
    mprintf("               to # atoms in target %i.\n",TargetMap.natom);
  }

  // Atoms have now been assigned IDs. Match up the unique strings in Ref with 
  // unique strings in target.
  AtomMap=(int*) malloc( RefMap.natom * sizeof(int));
  for (refatom=0; refatom<RefMap.natom; refatom++) {
    AtomMap[refatom]=-1;
    // If the ID of this reference atom is unique, look for same ID in target
    if (RefMap.M[refatom].isUnique) {
      for (targetatom=0; targetatom<TargetMap.natom; targetatom++) {
        // If ID of thie target atom is unique, check if it matches reference atom ID
        if (TargetMap.M[targetatom].isUnique) {
          if ( strcmp(TargetMap.M[targetatom].unique,
                      RefMap.M[refatom].unique)==0 ) {
            AtomMap[refatom]=targetatom;
            //mprintf("* TargetAtom %i(%s) maps to RefAtom %i(%s)\n",
            //        targetatom,TargetMap.P->names[targetatom],
            //        refatom,RefMap.P->names[refatom]);
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

  // DEBUG - Print current atom map and count # mapped atoms
  for (refatom=0; refatom<RefMap.natom; refatom++) {
    targetatom=AtomMap[refatom];
    if (targetatom>=0) {
      mprintf("* TargetAtom %6i(%4s) maps to RefAtom %6i(%4s)\n",
                      targetatom,TargetMap.P->names[targetatom],
                      refatom,RefMap.P->names[refatom]);
      numMappedAtoms++;
    } else {
      mprintf("* Could not map any TargetAtom to RefAtom %6i(%4s)\n",
                      refatom,RefMap.P->names[refatom]);
    }
  }
  mprintf("* %i total atoms were mapped.\n",numMappedAtoms);

  // Get list of atoms
  //RefMap.fourAtoms(fouratomlist);

  // Set up mapped atoms into frames for RMSD fitting
  refFrame = new Frame(numMappedAtoms,NULL);
  tgtFrame = new Frame(numMappedAtoms,NULL);
  refIndex = 0; // refIndex is atom number of the frames
  // DEBUG - targetIndex is refIndex*3
  targetIndex = 0;
  for (refatom=0; refatom<RefMap.natom; refatom++) {
    targetatom=AtomMap[refatom];
    if (targetatom>=0) {
      refFrame->SetCoord(refIndex, RefMap.F->Coord(refatom));
      tgtFrame->SetCoord(refIndex, TargetMap.F->Coord(targetatom));
      // DEBUG
      pdb_write_ATOM(buffer,PDBATOM,refIndex+1,TargetMap.P->names[targetatom],
                   (char*)"UNK",' ',1,tgtFrame->X[targetIndex],tgtFrame->X[targetIndex+1],
                   tgtFrame->X[targetIndex+2],0.0,0.0,(char*)"\0");
      targetIndex+=3;
      mprinterr("%s",buffer);
      // END DEBUG
      //mprintf("Ref %6i: ",refatom);
      //refFrame->printAtomCoord(refIndex);
      //mprintf("Tgt %6i: ",targetatom);
      //tgtFrame->printAtomCoord(refIndex);
      refIndex++;
    }
  }
  // Perform RMSD fit of target mapped atoms to reference mapped atoms
  rms = tgtFrame->RMSD(refFrame,Rot,Trans,false);
  mprintf("Rms of target map from ref map is %lf\n",rms);
  // Rotate and translate original target frame
  delete tgtFrame;
  tgtFrame = TargetMap.F->Copy();
  tgtFrame->Translate(Trans);
  tgtFrame->Rotate(Rot);
  tgtFrame->Translate(Trans+3);
  // DEBUG - write fit target out to PDB
  //targetIndex=0;
  //for (targetatom=0; targetatom<TargetMap.P->natom; targetatom++) {
  //  pdb_write_ATOM(buffer,PDBATOM,targetatom+1,TargetMap.P->names[targetatom],
  //                 (char*)"UNK",' ',1,tgtFrame->X[targetIndex],tgtFrame->X[targetIndex+1],
  //                 tgtFrame->X[targetIndex+2],0.0,0.0,(char*)"\0");
  //  targetIndex+=3;
  //  mprinterr("%s",buffer);
  //}

  // Cleanup
  free(AtomMap);
  delete refFrame;
  delete tgtFrame;

  return 0;
}

