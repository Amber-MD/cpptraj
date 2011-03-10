// AtomMap
#include <cstdlib>
#include <cstring>
#include "AtomMap.h"
#include "CpptrajStdio.h"
#include "TorsionRoutines.h"
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
        if ((nbondi<MAXBONDS)&&(nbondj<MAXBONDS)) {
          M[i].bond[nbondi++]=j;
          M[j].bond[nbondj++]=i;
          M[i].nbond=nbondi;
          M[j].nbond=nbondj;
          if (debug>1) mprintf("BONDED!\n");
        } else {
          mprintf("Warning: Bonding %s:%i to %s:%i; Valences MAXED (>%i)!\n",
                  names[i],i,names[j],j,MAXBONDS);
        }
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
 * atommap::markComplete()
 * Go through each atom in the map. If the atom is unique and all bonded
 * atoms are unique mark the atom as completely mapped.
 * Print bond information for each atom in the map, indicate whether
 * atoms are unique or not.
 */
void atommap::markComplete() {
  int atom,bond,bondatom,nunique;

  for (atom=0; atom<natom; atom++) {
    if (debug>0) mprintf("  Atom %4i: %s-%2i |",atom,names[atom],M[atom].isUnique);
    nunique=0;
    for (bond=0; bond<M[atom].nbond; bond++) {
      bondatom = M[atom].bond[bond];
      if (debug>0) mprintf(" %4i:%s-%2i",bondatom,names[bondatom],M[bondatom].isUnique);
      if (M[atom].isUnique && M[bondatom].isUnique) nunique++;
      if (nunique==M[atom].nbond) {
        if (debug>0) mprintf(" Atom is completely mapped.");
        M[atom].complete=true;
      }
    }
    
    if (debug>0) mprintf("\n");
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
  int atom,bond;
  char *ptr;

  natom=P->natom;
  names=(char**) malloc(natom * sizeof(char*));
  // Set up atom names.
  for (atom=0; atom<natom; atom++) {
    names[atom]=(char*) malloc(2*sizeof(char));
    // position ptr at first non-space character in name
    ptr=P->names[atom];
    while (*ptr==' ' && *ptr!='\0') ptr++;
    // if NULL something went wrong, abort
    if (*ptr=='\0') {
      strcpy(names[atom],"");
      continue;
    }
    names[atom][0]=ptr[0];
    // If C, check for L or l for chlorine
    if (ptr[0]=='C') {
      if (ptr[1]=='L' || ptr[1]=='l') names[atom][0]='X';
    }
    names[atom][1]='\0';
    // DEBUG
    mprintf("  Atom %i element: [%s]\n",atom,names[atom]);
  }
  // Allocate memory for atoms and initialize each atom
  M=(mapatom*) malloc( natom * sizeof(mapatom));
  for (atom=0; atom<natom; atom++) {
    for (bond=0; bond<MAXBONDS; bond++) M[atom].bond[bond]=-1;
    M[atom].nbond=0;
    M[atom].complete=false;
    memset(M[atom].atomID,' ',ATOMIDLENGTH);
    memset(M[atom].unique,' ',UNIQUELENGTH);
    M[atom].isUnique=1; // Assume unique until proven otherwise
  }
  return 0;
}

// ============================================================================
/*
 * AtomMap::mapChiral()
 * Given two atommaps and a map relating the two, find chiral centers for
 * which at least 3 of the atoms have been mapped. Assign the remaining
 * two atoms based on improper dihedrals. ONLY WORKS FOR SP3
 */
int AtomMap::mapChiral(atommap *Ref, atommap *Tgt, int *AMap) {
  int atom,tatom,bond,nunique,notunique;
  int uR[5], uT[5], r, t, nR[4], nT[4];
  double dR[4], dT[4], delta;
  bool newchiral;

  newchiral=true;

  while (newchiral) {

  // newchiral will be set to true only if new atoms get mapped.
  newchiral=false;

  for (atom=0; atom<Ref->natom; atom++) {
    if (!Ref->M[atom].isUnique) continue;
    tatom = AMap[atom];
    if (tatom<0) {
      mprintf("Error: AtomMap::mapChiral: Atom %i in reference is unique but not mapped!\n",
              atom);
      return 1;
    }
    // Sanity check - if Ref atom is completely mapped, target should be.
    // If so, skip
    if (Ref->M[atom].complete) {
      if (!Tgt->M[tatom].complete) {
        mprintf("Error: AtomMap::mapChiral: Ref atom %i is complete but Tgt atom %i is not.\n",
                atom,tatom);
        return 1;
      }
      continue;
    }
    // Check for SP3 (4 bonds)
    if (Ref->M[atom].nbond!=4) continue;
    if (Tgt->M[tatom].nbond!=4) {
      mprintf("Error: AtomMap:::mapChiral: Ref atom %i is SP3 but Tgt atom %i is not!\n",
              atom,tatom);
      return 1;
    }
    // Place bonded atoms (starting with central atom) in R and T
    uR[0] = atom;
    uT[0] = tatom;
    nunique=1;
    notunique=0;
    for (bond=0; bond<Ref->M[atom].nbond; bond++) {
      r = Ref->M[atom].bond[bond];
      t = AMap[r];
      if (Ref->M[r].isUnique && Tgt->M[t].isUnique) {
        uR[nunique] = r;
        uT[nunique] = t;
        nunique++;
      } else {
        nR[notunique] = r;
        nT[notunique] = t;
        notunique++;
      }
    }
    // At this point non-unique atoms in nT are -1 by definition. Fill nT
    // with nonunique atoms from target
    r=0;
    for (bond=0; bond<Tgt->M[tatom].nbond; bond++) {
      t = Tgt->M[tatom].bond[bond];
      if (!Tgt->M[t].isUnique) nT[r++] = t;
    }
    mprintf("  Potential Chiral center %i_%s/%i_%s: Unique atoms=%i, non-Unique=%i/%i\n",
            atom,Ref->names[atom],tatom,Tgt->names[tatom],nunique,notunique,r);
    for (r=0; r<nunique; r++)
      mprintf("\t   Mapped\t%4i %4i\n",uR[r],uT[r]);
    for (r=0; r<notunique; r++)
      mprintf("\tNotMapped\t%4i %4i\n",nR[r],nT[r]);
    // If all atoms are unique no need to map
    if (nunique==5) continue;
    // Require at least 3 unique atoms for dihedral calc. If < 3 unique this is
    // probably not really a chircal center anyway.
    if (nunique<3) {
      mprintf("    Warning: Center has < 3 mapped atoms, dihedral cannot be calcd.\n");
      continue;
    }
    // Calculate reference improper dihedrals
    for (r=0; r<notunique; r++) {
      dR[r] = Torsion(Ref->F->Coord(uR[0]),Ref->F->Coord(uR[1]),
                      Ref->F->Coord(uR[2]),Ref->F->Coord(nR[r]));
      mprintf("    Ref Improper %i = %lf\n",r,dR[r]);
      dT[r] = Torsion(Tgt->F->Coord(uT[0]),Tgt->F->Coord(uT[1]),
                      Tgt->F->Coord(uT[2]),Tgt->F->Coord(nT[r]));
      mprintf("    Tgt Improper %i = %lf\n",r,dT[r]);
    }
    // Match impropers to each other using a cutoff
    // NOTE: 10.0 seems reasonable? Also there is currently no check for 
    //       repeated deltas.
    for (r=0; r<notunique; r++) {
      for (t=0; t<notunique; t++) {
        delta = dR[r] - dT[t];
        if (delta<0.0) delta=-delta;
        if (delta<10.0) {
          mprintf("    Mapping tgt atom %i to ref atom %i based on chirality.\n",nR[r],nT[t]);
          AMap[ nR[r] ]=nT[t];
          // Once an atom has been mapped set its unique flag
          Ref->M[nR[r]].isUnique=1;
          Tgt->M[nT[t]].isUnique=1;
          // Since new atoms have been mapped, new chiral centers may have appeared
          newchiral=true;
        }
      }
    }
    // Check if atom is now completely mapped
    nunique=0;
    for (bond=0; bond < Ref->M[atom].nbond; bond++) {
      r = Ref->M[atom].bond[bond];
      t = AMap[r];
      if (t>-1) nunique++;
    }
    if (nunique==Ref->M[atom].nbond) {
      Ref->M[atom].complete=true;
      Tgt->M[tatom].complete=true;
    } 
  } // End loop over natom

  } // End loop over newchiral

  return 0;
}

/*
 * AtomMap::mapByIndex()
 * Given to atommaps and a map relating the two, attempt to map any remaining
 * incomplete atoms by assuming the atom indices in reference and target are
 * in similar orders. At this point all unique atoms and chiral centers should
 * have been identified, and the remaining unmapped atoms are those in very
 * symmetrical regions like methyl groups, rings etc.
 */
int AtomMap::mapByIndex(atommap *Ref, atommap *Tgt, int *AMap) {
  int atom,tatom,bond,tbond,numMapped,iterations,r,t,nunique;
  bool allAtomsMapped=false;

  iterations=0;
  while (!allAtomsMapped) {
    numMapped=0;
    // Find next incompletely mapped atom
    for (atom=0; atom<Ref->natom; atom++) {
      // Skip over non-unique (non-mapped) atoms
      if (!Ref->M[atom].isUnique) continue;
      tatom = AMap[atom];
      if (tatom<0) {
        mprintf("Error: AtomMap::mapByIndex: Atom %i in reference is unique but not mapped!\n",
                atom);
        return 1;
      }
      // Check that num bonds match in Ref and target
      if (Ref->M[atom].nbond!=Tgt->M[tatom].nbond) {
        mprintf("Error: AtomMap:::mapByIndex: Ref atom %i #bonds does not match Tgt atom %i\n",
                atom,tatom);
        return 1;
      }
      // Skip completely mapped atoms - check that both Ref and Tgt are complete
      if (Ref->M[atom].complete) {
        if (!Tgt->M[tatom].complete) {
          mprintf("Error: AtomMap:mapByIndex: Ref atom %i is complete but Tgt atom %i is not.\n",
                atom,tatom);
          return 1;
        }
        numMapped++;
        continue;
      }
      // This atom is mapped, but bonded atoms are not completely mapped. Map 
      // remaining unmapped atoms in target to unmapped atoms in reference. 
      // Map first by name, then by index.
      for (bond=0; bond < Ref->M[atom].nbond; bond++) {
        r = Ref->M[atom].bond[bond];
        if (Ref->M[r].isUnique) continue;
        for (tbond=0; tbond < Tgt->M[tatom].nbond; tbond++) {
          t = Tgt->M[tatom].bond[tbond];
          if (Tgt->M[t].isUnique) continue;
          // Atom r bonded to atom, atom t bonded to tatom. r and t are not
          // yet mapped. Check if names match
          if (strcmp(Ref->names[r],Tgt->names[t])!=0) continue;
          // Here we could check to make sure the index order matches as well
          // MAP THEM
          mprintf("  Mapping Tgt %i to Ref %i based on name/bonding.\n",t,r);
          AMap[r] = t;
          Ref->M[r].isUnique=1;
          Tgt->M[t].isUnique=1;
          // Since r has been mapped, exit the loop to avoid remapping
          break;
        }
      }
      // Check if atom is completely mapped now
      nunique=0;
      for (bond=0; bond < Ref->M[atom].nbond; bond++) {
        r = Ref->M[atom].bond[bond];
        t = AMap[r];
        if (t>-1) nunique++;
      }
      if (nunique==Ref->M[atom].nbond) {
        Ref->M[atom].complete=true;
        Tgt->M[tatom].complete=true;
        numMapped++;
      }
    } // End loop over atoms
    if (numMapped >= Ref->natom) allAtomsMapped=true;
    // Safety valve - shouldnt have to map more times than # atoms
    if (iterations > Ref->natom) {
      mprintf("Error: AtomMap::mapByIndex: # iterations greater than # atoms.\n");
      return 1;
    }
    iterations++;
  } // End - allAtomsMapped
  return 0;
}

/*
 * AtomMap::init()
 * Expected call: atommap <target> <reference>
 * Attempt to create a map from atoms in target to atoms in reference solely
 * based on how they are bonded (not how they are named). 
 */
int AtomMap::init() {
  char *refName, *targetName;
  int refIndex, targetIndex;
  int refatom,targetatom, *AMap;
  int numMappedAtoms=0;
  // For RMSD
  Frame *refFrame, *tgtFrame;
  double rms, Rot[9],Trans[6];
  // DEBUG - pdb write
  char buffer[83];
  PtrajFile fitpdb;
  PtrajFile mappdb;
  // END DEBUG

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

  TargetMap.setup();
  TargetMap.calcDist();
  TargetMap.determineAtomID();

  // Number of atoms in each map MUST be equal
  if (RefMap.natom!=TargetMap.natom) {
    mprintf("AtomMap::init: Warning: # atoms in reference %i not equal\n",RefMap.natom);
    mprintf("               to # atoms in target %i.\n",TargetMap.natom);
  }

  // Atoms have now been assigned IDs. Match up the unique strings in Ref with 
  // unique strings in target.
  AMap=(int*) malloc( RefMap.natom * sizeof(int));
  for (refatom=0; refatom<RefMap.natom; refatom++) {
    AMap[refatom]=-1;
    // If the ID of this reference atom is unique, look for same ID in target
    if (RefMap.M[refatom].isUnique) {
      for (targetatom=0; targetatom<TargetMap.natom; targetatom++) {
        // If ID of thie target atom is unique, check if it matches reference atom ID
        if (TargetMap.M[targetatom].isUnique) {
          if ( strcmp(TargetMap.M[targetatom].unique,
                      RefMap.M[refatom].unique)==0 ) {
            // Check that number of bonds is consistent
            if (RefMap.M[refatom].nbond!=TargetMap.M[targetatom].nbond) {
              mprintf("Warning: AtomMap: Atoms R%i and T%i have same ID but different # bonds!\n",
                      refatom,targetatom);
            }
            AMap[refatom]=targetatom;
            //mprintf("* TargetAtom %i(%s) maps to RefAtom %i(%s)\n",
            //        targetatom,TargetMap.P->names[targetatom],
            //        refatom,RefMap.P->names[refatom]);
          } // If unique strings match
        } // If target atom is unique
      } // Loop over target atoms
    } // If reference atom is unique
  } // Loop over reference atoms

  // Search for completely mapped atoms. If an atom and all atoms
  // it is bonded to are unique, mark the atom as completely mapped.
  RefMap.markComplete();
  TargetMap.markComplete();

  // Map any unmapped chiral centers
  mapChiral(&RefMap, &TargetMap, AMap);
  //RefMap.markComplete();
  //TargetMap.markComplete();

  // Try to map the rest by atom name and connectivity
  mapByIndex(&RefMap, &TargetMap, AMap);

  // DEBUG - Print current atom map and count # mapped atoms
  for (refatom=0; refatom<RefMap.natom; refatom++) {
    targetatom=AMap[refatom];
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
  if (numMappedAtoms!=RefMap.natom)
    mprintf("Warning: Not all atoms were mapped.\n");

  // Get list of atoms
  //RefMap.fourAtoms(fouratomlist);

  // Set up mapped atoms into frames for RMSD fitting
  refFrame = new Frame(numMappedAtoms,NULL);
  tgtFrame = new Frame(numMappedAtoms,NULL);
  refIndex = 0; // refIndex is atom number of the frames
  // DEBUG - For PDB write, targetIndex is refIndex*3
  targetIndex = 0;
  mappdb.SetupFile((char*)"targetmap.pdb\0",WRITE,PDBFILE,STANDARD,0);
  mappdb.OpenFile();
  fitpdb.SetupFile((char*)"fit.pdb\0",WRITE,PDBFILE,STANDARD,0);
  fitpdb.OpenFile();
  // END DEBUG
  for (refatom=0; refatom<RefMap.natom; refatom++) {
    targetatom=AMap[refatom];
    if (targetatom>=0) {
      refFrame->SetCoord(refIndex, RefMap.F->Coord(refatom));
      tgtFrame->SetCoord(refIndex, TargetMap.F->Coord(targetatom));
      // DEBUG
      pdb_write_ATOM(buffer,PDBATOM,refIndex+1,TargetMap.P->names[targetatom],
                   (char*)"UNK",' ',1,tgtFrame->X[targetIndex],tgtFrame->X[targetIndex+1],
                   tgtFrame->X[targetIndex+2],0.0,0.0,(char*)"\0");
      targetIndex+=3;
      mappdb.IO->Printf("%s",buffer);
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
  targetIndex=0;
  for (targetatom=0; targetatom<TargetMap.P->natom; targetatom++) {
    pdb_write_ATOM(buffer,PDBATOM,targetatom+1,TargetMap.P->names[targetatom],
                   (char*)"UNK",' ',1,tgtFrame->X[targetIndex],tgtFrame->X[targetIndex+1],
                   tgtFrame->X[targetIndex+2],0.0,0.0,(char*)"\0");
    targetIndex+=3;
    fitpdb.IO->Printf("%s",buffer);
  }
  mappdb.CloseFile();
  fitpdb.CloseFile();
  // END DEBUG

  // Cleanup
  free(AMap);
  delete refFrame;
  delete tgtFrame;

  return 0;
}

