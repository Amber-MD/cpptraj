// AtomMap
#include <cstdlib>
#include <cstring>
#include "AtomMap.h"
#include "CpptrajStdio.h"
#include "TorsionRoutines.h"
// DEBUG
//#include "Mol2File.h"

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

// atommap CONSTRUCTOR
atommap::atommap() {
  M=NULL;
  natom=0;
  names=NULL;
  F=NULL;
  debug=0;
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

// Set atommap debug
void atommap::SetDebug(int debugIn) {
  debug=debugIn;
}

/*
 * atommap::getCut() 
 * Return a cutoff based on optimal covalent bond distance based on the 
 * identities of atom1 and atom2. When multiple hybridizations are possible
 * the longest possible bond length is used.
 * Treat X as chlorine for now, Y as Bromine.
 * Unless otherwise noted values taken from:
 *   Huheey, pps. A-21 to A-34; T.L. Cottrell, "The Strengths of Chemical Bonds," 
 *       2nd ed., Butterworths, London, 1958; 
 *   B. deB. Darwent, "National Standard Reference Data Series," National Bureau of Standards, 
 *       No. 31, Washington, DC, 1970; S.W. Benson, J. Chem. Educ., 42, 502 (1965).
 * Can be found on the web at:
 *   http://www.wiredchemist.com/chemistry/data/bond_energies_lengths.html
 */
double atommap::getCut(char *atom1, char *atom2) {
  // Default cutoff
  double cut=1.60;

  // Self
  if (strcmp(atom1,atom2)==0) {
    if      (strcmp(atom1,"H")==0) cut=0.74;
    else if (strcmp(atom1,"N")==0) cut=1.45;
    else if (strcmp(atom1,"C")==0) cut=1.54;
    else if (strcmp(atom1,"O")==0) cut=1.48;
    else if (strcmp(atom1,"P")==0) cut=2.21;
    else if (strcmp(atom1,"S")==0) cut=2.05; // Gas-phase value, S=S is 1.49
  }
  // Bonds to H 
  else if ( compareName(atom1,atom2,"H","C")==0 )
    cut=1.09;
  else if ( compareName(atom1,atom2,"H","N")==0 )
    cut=1.01;
  else if ( compareName(atom1,atom2,"H","O")==0 )
    cut=0.96;
  else if ( compareName(atom1,atom2,"H","P")==0 )
    cut=1.44;
  else if ( compareName(atom1,atom2,"H","S")==0 )
    cut=1.34;
  // Bonds to C
  else if ( compareName(atom1,atom2,"C","N")==0 )
    cut=1.47;
  else if ( compareName(atom1,atom2,"C","O")==0 )
    cut=1.43;
  else if ( compareName(atom1,atom2,"C","P")==0 )
    cut=1.84;
  else if ( compareName(atom1,atom2,"C","F")==0 )
    cut=1.35;
  else if ( compareName(atom1,atom2,"C","X")==0 )
    cut=1.77;
  else if ( compareName(atom1,atom2,"C","Y")==0 )
    cut=1.94;
  else if ( compareName(atom1,atom2,"C","S")==0 )
    cut=1.82;
  // Bonds to N
  else if ( compareName(atom1,atom2,"N","O")==0 )
    cut=1.40;
  else if ( compareName(atom1,atom2,"N","S")==0 )
    cut=1.68; // Postma & Vos, Acta Cryst. (1973) B29, 915
  else if ( compareName(atom1,atom2,"N","F")==0 )
    cut=1.36;
  else if ( compareName(atom1,atom2,"N","X")==0 )
    cut=1.75;

  // Bonds to P
  else if ( compareName(atom1,atom2,"P","O")==0 )
    cut=1.63;
  else if ( compareName(atom1,atom2,"P","S")==0 )
    cut=1.86;
  else if ( compareName(atom1,atom2,"P","F")==0 )
    cut=1.54;
  else if ( compareName(atom1,atom2,"P","X")==0 )
    cut=2.03; 

  // Bonds to O
  else if ( compareName(atom1,atom2,"O","S")==0 )
    cut=1.48;
  else if ( compareName(atom1,atom2,"O","F")==0 )
    cut=1.42;

  // Bonds to S
  else if ( compareName(atom1,atom2,"S","F")==0 )
    cut=1.56;
  else if ( compareName(atom1,atom2,"S","X")==0 )
    cut=2.07;

  // No cutoff, use default
  else {
    if (debug>0) {
      mprintf("Warning: atommap::getCut: Cut not found for %s - %s\n",atom1,atom2);
      mprintf("                          Using default cutoff of %lf\n",cut);
    }
  }

  // Padding value
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
  // DEBUG
  //bool isRepeated;
  //int k, atom2;

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

  // DEBUG
  // For each atom with a truly unique ID, determine if it is bonded to a
  // non-unique partner. If that partner is itself unique among bonded
  // partners (e.g. H2-C-N where C is unique, N is unique by extension),
  // give it a unique ID of atomID-element
/*  for (i=0; i<natom; i++) {
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
            mprintf("DBG: Non-unique Atom %i:%s could be unique by extension.\n",
                    atom,names[atom]);
            //strcpy(M[atom].unique,M[i].unique);
            //strcat(M[atom].unique,"-");
            //strcat(M[atom].unique,names[atom]);
            //M[atom].isUnique=1;
          }
        } // End if bonded atom is not unique
      } // End loop j over bonds of atom i
    } // End if atom i is unique
  } // End loop i over atoms in map 
*/

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
 * atommap::BondIsRepeated()
 * Check if the name of the atom in bond <bond> bonded to atom <atom>
 * is the same as the name of any other non-unique atom bonded to 
 * atom <atom>.
 */
bool atommap::BondIsRepeated(int atom, int bond) {
  int bondedAtom,bondedAtom2;
  // If 1 or no bonds, atom cant possibly be repeated.
  if (M[atom].nbond<2) return false;
  bondedAtom = M[atom].bond[bond];
  for (int n=0; n < M[atom].nbond; n++) {
    if (n==bond) continue;
    bondedAtom2 = M[atom].bond[n];
    //mprintf("              %i) %i:%s %i:%s\n",n,bondedAtom,names[bondedAtom],
    //        bondedAtom2,names[bondedAtom2]);
    // If bondedAtom2 is unique dont check it
    if (M[bondedAtom2].isUnique) continue;
    if (strcmp(names[bondedAtom],names[bondedAtom2])==0) return true;
  }
  return false;
}

/* 
 * atommap::setup()
 * Allocate memory for atom map. In order to easily create the uniqueID 
 * strings the atom names need to be 1 char long. Convert chlorine 
 * to X for now, bromine to Y etc.
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
    // If B, check for R or r for bromine
    if (ptr[0]=='B') {
      if (ptr[1]=='R' || ptr[1]=='r') names[atom][0]='Y';
    }
    names[atom][1]='\0';
    // DEBUG
    if (debug>0) mprintf("  Atom %i element: [%s]\n",atom,names[atom]);
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

// DEBUG
/*
void atommap::WriteMol2(char *m2filename) {
  Mol2File outfile;
  int *Selected;
  // Temporary parm to play with
  AmberParm *tmpParm;

  // Create mask containing all atoms
  Selected = (int*) malloc(natom*sizeof(int));
  for (int atom=0; atom<natom; atom++) Selected[atom]=atom;
  // Fake strip, just use as crap way to copy
  tmpParm = P->modifyStateByMask(Selected,natom);
  free(Selected);
  // Modify the bonds array to include this info
  tmpParm->ResetBondInfo();
  for (int atom=0; atom<natom; atom++) 
    for (int bond=0; bond < M[atom].nbond; bond++) 
      tmpParm->AddBond(atom, M[atom].bond[bond], 0);

  // Trajectory Setup
  outfile.File=new PtrajFile();
  outfile.File->SetupFile(m2filename,WRITE,MOL2FILE,UNKNOWN_TYPE,debug);
  outfile.trajfilename = outfile.File->basefilename;
  outfile.debug=debug;
  outfile.SetTitle(m2filename);
  outfile.P=P;
  outfile.SetupWrite();
  outfile.open();
  outfile.F=F;
  outfile.writeFrame(0);
  outfile.close();
  delete tmpParm;
}*/
// ============================================================================

// CONSTRUCTOR
AtomMap::AtomMap() {
  AMap=NULL;
  newFrame=NULL;
  newParm=NULL;
  stripParm=NULL;
  maponly=false;
}

// DESTRUCTOR
AtomMap::~AtomMap() {
  if (AMap!=NULL) free(AMap);
  if (newFrame!=NULL) delete newFrame;
  if (newParm!=NULL) delete newParm;
  if (stripParm!=NULL) delete stripParm;
}

/*
 * AtomMap::mapSingleBonds()
 */
int AtomMap::mapSingleBonds(atommap *Ref, atommap *Tgt) {
  int atom,bond,r;
  int tatom,tbond,t;

  for (atom=0; atom < Ref->natom; atom++) {
    // Skip non-unique atoms in Ref
    if (!Ref->M[atom].isUnique) continue;
    // Unique atoms should be mapped
    tatom = AMap[atom];
    if (tatom<0) {
      mprintf("      Warning: mapSingleBonds: Ref %i is unique but not mapped.\n",atom);
      continue;
    }
    // For each non-unique atom bonded to unique Reference atom, try to 
    // find a matching non-unique atom in Target by virtue of it being the
    // only possible match.
    for (bond=0; bond < Ref->M[atom].nbond; bond++) {
      r = Ref->M[atom].bond[bond];
      // Check that bonded atom r is not unique
      if (Ref->M[r].isUnique) continue;
      //mprintf("        Ref: Checking non-unique %i:%s bonded to %i:%s\n",r,Ref->names[r],
      //        atom,Ref->names[atom]);
      // Check that non-unique bonded ref atom r name is not the same as any other 
      // non-unique bonded atom.
      if ( Ref->BondIsRepeated(atom, bond) ) continue;
      // At this point r is the only one of its kind bonded to atom.
      // Check if there is an analogous atom bonded to unique Target atom
      // tatom.
      for (tbond=0; tbond < Tgt->M[tatom].nbond; tbond++) {
        t = Tgt->M[tatom].bond[tbond];
        // Check that bonded atom t is not unique
        if (Tgt->M[t].isUnique) continue;
        //mprintf("          Tgt: Checking non-unique %i:%s bonded to %i:%s\n",t,Tgt->names[t],
        //        tatom,Tgt->names[tatom]);
        // Check that non-unique bonded tgt atom t name is not the same as any other
        // non-unique bonded atom.
        if ( Tgt->BondIsRepeated(tatom, tbond) ) continue;
        // At this point t is the only one of its kind bonded to tatom.
        // Check if its name matches r. If so, map it.
        if ( strcmp(Ref->names[r], Tgt->names[t])==0 ) {
          if (debug>0) 
            mprintf("    Mapping tgt %i to ref %i based on single bond to unique.\n",t,r);
          AMap[r]=t;
          Ref->M[r].isUnique=1;
          Tgt->M[t].isUnique=1;
        }
      } // End loop over atoms bonded to tatom
    } // End loop over atoms bonded to atom
  } // End loop over ref atoms
  return 0;
}        

/*
 * AtomMap::mapChiral()
 * Given two atommaps and a map relating the two, find chiral centers for
 * which at least 3 of the atoms have been mapped. Assign the remaining
 * two atoms based on improper dihedrals. ONLY WORKS FOR SP3
 */
int AtomMap::mapChiral(atommap *Ref, atommap *Tgt) {
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
      mprintf("      Warning: mapChiral: Atom %i in reference is unique but not mapped!\n",
              atom);
      continue;
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
      // Bonded atom r is not unique
      if (!Ref->M[r].isUnique) {
        nR[notunique++] = r;
      // Bonded atom r is unique. If a target was mapped to it and is unique
      // (i.e. it is the same atom) store it.
      } else {
        if (t>=0) {
          if (Ref->M[r].isUnique && Tgt->M[t].isUnique) {
            uR[nunique] = r;
            uT[nunique] = t;
            nunique++;
          }
        } 
      }
    }
    // At this point non-unique atoms in nT are -1 by definition. Fill nT
    // with nonunique atoms from target
    r=0;
    for (bond=0; bond<Tgt->M[tatom].nbond; bond++) {
      t = Tgt->M[tatom].bond[bond];
      if (!Tgt->M[t].isUnique) nT[r++] = t;
    }
    // notunique should be the same as r
    if (notunique!=r) mprintf("Warning: Ref and Tgt do not have the same # of nonunique atoms.\n");
    if (debug>0) { 
      mprintf("  Potential Chiral center %i_%s/%i_%s: Unique atoms=%i, non-Unique=%i/%i\n",
              atom,Ref->names[atom],tatom,Tgt->names[tatom],nunique,notunique,r);
      for (r=0; r<nunique; r++)
        mprintf("\t   Mapped\t%4i %4i\n",uR[r],uT[r]);
      for (r=0; r<notunique; r++)
        mprintf("\tNotMapped\t%4i %4i\n",nR[r],nT[r]);
    }
    // If all atoms are unique no need to map
    if (nunique==5) continue;
    // Require at least 3 unique atoms for dihedral calc. If < 3 unique this is
    // probably not really a chircal center anyway.
    if (nunique<3) {
      if (debug>0) 
        mprintf("    Warning: Center has < 3 mapped atoms, dihedral cannot be calcd.\n");
      continue;
    }
    // Calculate reference improper dihedrals
    for (r=0; r<notunique; r++) {
      dR[r] = Torsion(Ref->F->Coord(uR[0]),Ref->F->Coord(uR[1]),
                      Ref->F->Coord(uR[2]),Ref->F->Coord(nR[r]));
      if (debug>1) mprintf("    Ref Improper %i = %lf\n",r,dR[r]);
      dT[r] = Torsion(Tgt->F->Coord(uT[0]),Tgt->F->Coord(uT[1]),
                      Tgt->F->Coord(uT[2]),Tgt->F->Coord(nT[r]));
      if (debug>1) mprintf("    Tgt Improper %i = %lf\n",r,dT[r]);
    }
    // Match impropers to each other using a cutoff
    // NOTE: 10.0 seems reasonable? Also there is currently no check for 
    //       repeated deltas.
    for (r=0; r<notunique; r++) {
      for (t=0; t<notunique; t++) {
        delta = dR[r] - dT[t];
        if (delta<0.0) delta=-delta;
        if (delta<10.0) {
          if (debug>0)
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
 * AtomMap::mapUniqueRefToTgt()
 * If the number of atoms in Ref is different from Tgt, it is possible that
 * Tgt is missing atoms (or maybe vice versa). If the difference is not too
 * great it may be possible to look for an unmapped atom in Ref that has
 * same name and at least 1 matching bond (# bonds may be diff due to # atoms).
 * If only one name matches, probably safe to map it. Return 1 if the atom
 * could be mapped, 0 otherwise.
 */
int AtomMap::mapUniqueRefToTgt(atommap *Ref, atommap *Tgt, int atom) {
  int t,commonBond,bond,tbond,match,r;
  bool alreadyMapped;

  match=-1;
  for (t=0; t < Tgt->natom; t++) {
    //mprintf("DBG:        %i:%s %i:%s\n",t,Tgt->names[t],atom,Ref->names[atom]);
    // If atom #s are different Tgt atom could be unique but not mapped. Check
    // if tgt has already been mapped using Amap. 
    alreadyMapped=false;
    for (bond=0; bond < Ref->natom; bond++) {
      if (AMap[bond]==t) {alreadyMapped=true; break;}
    }
    if (alreadyMapped) continue;
    // Check name
    if ( strcmp(Tgt->names[t], Ref->names[atom])==0) {
      if (debug>1) mprintf("        Attempting match of %i:%s to %i:%s\n",
              t,Tgt->P->names[t],atom,Ref->P->names[atom]);
      // Check that at least 1 bond is in common
      commonBond=0;
      for (bond=0; bond < Ref->M[atom].nbond; bond++) {
        // Check Map for ref bonded atom
        r = AMap[ Ref->M[atom].bond[bond] ];
        if (debug>1) 
          mprintf("          Ref %i:%s bonded to %i:%s (%i in tgt)\n",atom,Ref->P->names[atom],
                  Ref->M[atom].bond[bond], Ref->P->names[ Ref->M[atom].bond[bond] ],r);
        // If no mapping exists cant check it
        if (r<0) continue;
        for (tbond=0; tbond < Tgt->M[t].nbond; tbond++) {
          if (debug>1)
            mprintf("            Tgt %i:%s bonded to %i:%s\n",t,Tgt->P->names[t],
                    Tgt->M[t].bond[tbond], Tgt->P->names[ Tgt->M[t].bond[tbond] ]);
          if (r == Tgt->M[t].bond[tbond]) 
            commonBond++;
        }
      }
      if (commonBond==0) continue;
      // This Tgt Name matches and at least 1 bond in common with Ref atom
      // Check that a match has not yet been found for ref
      if (match!=-1) {
        mprintf("      Warning: mapUniqueRefToTgt: Ref %i:%s has multiple potential matches\n",
                atom,Ref->P->names[atom]);
        mprintf("               among Tgt [%i:%s, %i:%s]\n",
                t,Tgt->P->names[t],match,Tgt->P->names[match]);
        return 0;
      }
      match = t;
    }
  }
  if (match==-1) return 0;
  // Only one match found - map it
  if (debug>0) 
    mprintf("    Mapping target %i:%s to unique ref %i:%s\n",match,Tgt->P->names[match],
            atom,Ref->P->names[atom]);
  AMap[atom]=match;
  Tgt->M[match].isUnique=1;
  return 1;
}

/*
 * AtomMap::mapByIndex()
 * Given to atommaps and a map relating the two, attempt to map any remaining
 * incomplete atoms by assuming the atom indices in reference and target are
 * in similar orders. At this point all unique atoms and chiral centers should
 * have been identified, and the remaining unmapped atoms are those in very
 * symmetrical regions like methyl groups, rings etc.
 */
int AtomMap::mapByIndex(atommap *Ref, atommap *Tgt) {
  int atom,tatom,bond,tbond,numMapped,iterations,r,t,nunique,lastMapped;
  bool allAtomsMapped=false;

  iterations=0;
  lastMapped=-1;
  while (!allAtomsMapped) {
    numMapped=0;
    // Find next incompletely mapped atom
    for (atom=0; atom<Ref->natom; atom++) {
      // Skip over non-unique (non-mapped) atoms
      if (!Ref->M[atom].isUnique) continue;
      tatom = AMap[atom];
      if (tatom<0) {
        mprintf("      Warning: mapByIndex: Atom %i in reference is unique but not mapped!\n",
                atom);
        // Reference atom is unique, but hasnt had a target mapped to it.
        // This can arise when the number of atoms in ref and tgt not equal. 
        // If the difference in atoms is not too great (probably 1) attempt
        // to look for a similar atom in tgt (name and index) and map it.
        if (mapUniqueRefToTgt(Ref,Tgt,atom)) numMapped++;
        continue;
      }
      // Check that num bonds match in Ref and target.
      // The # of bonds might not be equal if the # atoms in ref and tgt
      // not equal.
      if (Ref->M[atom].nbond!=Tgt->M[tatom].nbond) {
        mprintf("      Warning: mapByIndex: Ref atom %i #bonds %i does not match Tgt atom %i (%i)\n",
                atom,Ref->M[atom].nbond,tatom,Tgt->M[tatom].nbond);
        //return 1;
      }
      // Skip completely mapped atoms - check that both Ref and Tgt are complete
      if (Ref->M[atom].complete) {
        if (!Tgt->M[tatom].complete) {
          mprintf("Error: AtomMap::mapByIndex: Ref atom %i is complete but Tgt atom %i is not.\n",
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
          if (debug>0) mprintf("    Mapping Tgt %i to Ref %i based on name/bonding.\n",t,r);
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
    if (debug>1)
      mprintf("  mapByIndex: iteration %i: %i atoms mapped.\n",iterations,numMapped);
    // If we didnt map any more atoms this iteration than last, exit
    if (lastMapped == numMapped) break;
    lastMapped = numMapped;
    iterations++;
  } // End - allAtomsMapped
  return 0;
}

/*
 * AtomMap::init()
 * Expected call: atommap <target> <reference> [mapout <filename>] [maponly]
 * Attempt to create a map from atoms in target to atoms in reference solely
 * based on how they are bonded (not how they are named). 
 */
int AtomMap::init() {
  char *refName, *targetName, *outputname;
  PtrajFile outputfile;
  int refIndex, targetIndex;
  int refatom,targetatom;
  int numMappedAtoms=0;
  AtomMask *M1;
  
  RefMap.SetDebug(debug);
  TargetMap.SetDebug(debug);

  // Get Args
  outputname=A->getKeyString("mapout",NULL);
  maponly = A->hasKey("maponly");

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

  mprintf("    ATOMMAP: Atoms in trajectories associated with parm %s will be\n",
          TargetMap.P->parmName);
  mprintf("             mapped according to parm %s.\n",RefMap.P->parmName);
  if (outputname!=NULL)
    mprintf("             Map will be written to %s\n",outputname);
  if (maponly)
    mprintf("             maponly: Map will only be written, not used in trajectory read.\n");

  // For each map, set up (get element for each atom, initialize map mem),
  // determine what atoms are bonded to each other via simple distance
  // cutoffs, the give each atom an ID based on what atoms are bonded to
  // it, noting which IDs are unique for that map. 

  RefMap.setup();
  RefMap.calcDist();
  //RefMap.WriteMol2((char*)"RefMap.mol2\0"); // DEBUG
  RefMap.determineAtomID();

  TargetMap.setup();
  TargetMap.calcDist();
  //TargetMap.WriteMol2((char*)"TargetMap.mol2\0"); // DEBUG
  TargetMap.determineAtomID();

  // Number of atoms in each map MUST be equal
  if (RefMap.natom!=TargetMap.natom) {
    mprintf("      AtomMap::init: Warning: # atoms in reference (%i) not equal\n",RefMap.natom);
    mprintf("                     to # atoms in target (%i).\n",TargetMap.natom);
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
              mprintf("      Warning: AtomMap: Atoms R%i and T%i have same ID but different # bonds!\n",
                      refatom,targetatom);
            }
            AMap[refatom]=targetatom;
            if (debug>0) {
              mprintf("    Mapping Tgt %i(%s) to Ref %i(%s)\n",
                      targetatom,TargetMap.P->names[targetatom],
                      refatom,RefMap.P->names[refatom]);
            }
          } // If unique strings match
        } // If target atom is unique
      } // Loop over target atoms
    } // If reference atom is unique
  } // Loop over reference atoms

  // map atoms that are unique by virtue of being the only ones bonded to
  // other unique atoms
  mapSingleBonds(&RefMap, &TargetMap);

  // Search for completely mapped atoms. If an atom and all atoms
  // it is bonded to are unique, mark the atom as completely mapped.
  RefMap.markComplete();
  TargetMap.markComplete();

  // Map any unmapped chiral centers
  mapChiral(&RefMap, &TargetMap);
  //RefMap.markComplete();
  //TargetMap.markComplete();

  // Try to map the rest by atom name and connectivity
  mapByIndex(&RefMap, &TargetMap);

  // Print atom map and count # mapped atoms
  outputfile.SetupFile(outputname,WRITE,DATAFILE,UNKNOWN_TYPE,debug);
  outputfile.OpenFile();
  outputfile.IO->Printf("%-6s %4s %6s %4s\n","#TgtAt","Tgt","RefAt","Ref");
  for (refatom=0; refatom<RefMap.natom; refatom++) {
    targetatom=AMap[refatom];
    if (targetatom < 0) 
      outputfile.IO->Printf("%6s %4s %6i %4s\n","---","---",refatom+1,RefMap.P->names[refatom]);
    else
      outputfile.IO->Printf("%6i %4s %6i %4s\n",targetatom+1,TargetMap.P->names[targetatom],
                            refatom+1,RefMap.P->names[refatom]);
    if (targetatom>=0) {
      //mprintf("* TargetAtom %6i(%4s) maps to RefAtom %6i(%4s)\n",
      //                targetatom,TargetMap.P->names[targetatom],
      //                refatom,RefMap.P->names[refatom]);
      numMappedAtoms++;
    } //else {
    //  mprintf("* Could not map any TargetAtom to RefAtom %6i(%4s)\n",
    //                  refatom,RefMap.P->names[refatom]);
    //}
  }
  outputfile.CloseFile();
  mprintf("      %i total atoms were mapped.\n",numMappedAtoms);
  // If not all atoms could be mapped, return error since
  // this would screw up frame modification.
  // NOTE: Could do some sort of strip to remove atoms that
  // could not be mapped.
  if (numMappedAtoms!=RefMap.natom) {
    // If the number of mapped atoms is less than the number of reference
    // atoms but equal to the number of target atoms, can modify the reference
    // frame to only include mapped atoms
    if (numMappedAtoms<RefMap.natom && numMappedAtoms==TargetMap.natom) {
      // Create mask that includes only reference atoms that could be mapped
      M1 = new AtomMask();
      for (refatom=0; refatom<RefMap.natom; refatom++) {
        if (AMap[refatom]!=-1) M1->AddAtom(refatom);
      }
      // Strip reference parm
      mprintf("    Modifying reference %s topology and frame to match mapped atoms.\n",
              FL->FrameName(refIndex));
      stripParm = RefMap.P->modifyStateByMask(M1->Selected, numMappedAtoms);
      // Strip reference frame
      newFrame = new Frame(numMappedAtoms,RefMap.P->mass);
      newFrame->SetFrameFromMask(RefMap.F, M1);
      delete M1;
      // Replace reference with stripped versions
      if (FL->Replace(refIndex, newFrame, stripParm)) {
        mprintf("Error: AtomMap: Could not strip reference.\n");
        return 1;
      }
    } else {
      mprintf("Error: AtomMap: Not all atoms were mapped.\n");
      return 1;
    }
  }

  if (!maponly) {
    // Set up new Frame
    newFrame = new Frame(TargetMap.natom,TargetMap.P->mass);

    // Set up new Parm
    newParm = TargetMap.P->modifyStateByMap(AMap);
  }

  return 0;
}

/*
 * AtomMap::setup()
 * If the current parm does not match the target parm, deactivate. Otherwise
 * replace current parm with mapped parm.
 */
int AtomMap::setup() {
  if (maponly) {
    mprintf("    ATOMMAP: maponly was specified, not using atom map during traj read.\n");
    return 0;
  }
  if (P->pindex!=TargetMap.P->pindex ||
      P->natom !=TargetMap.P->natom) {
    if (debug>0) {
      mprintf("    ATOMMAP: Map for parm %s -> %s (%i atom).\n",TargetMap.P->parmName,
              RefMap.P->parmName,TargetMap.P->natom);
      mprintf("             Current parm %s (%i atom).\n",P->parmName,P->natom);
      mprintf("             Not using map for this parm.\n");
    }
    return 1;
  }
  mprintf("    ATOMMAP: Map for parm %s -> %s (%i atom).\n",TargetMap.P->parmName,
          RefMap.P->parmName,TargetMap.P->natom);

  P = newParm;
  
  return 0;
}

/*
 * AtomMap::action()
 * Modify the current frame based on the atom map. 
 */
int AtomMap::action() {
  if (maponly) return 0;
  for (int atom=0; atom < F->natom; atom++) 
    newFrame->SetCoord(atom, F->Coord(AMap[atom]));
  F = newFrame;
  return 0;
}

