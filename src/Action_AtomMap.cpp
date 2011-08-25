// AtomMap
#include <cstdlib>
#include <cstring>
#include "Action_AtomMap.h"
#include "CpptrajStdio.h"
#include "TorsionRoutines.h"
// DEBUG
#include "TrajectoryFile.h"
#include <cstdio>

//--------- PRIVATE ROUTINES ---------------------------------------
/* compareName(A1,A2,B1,B2) 
 * Compares pairs of names, return 1 if they match.  e.g. (A,B)==(A,B)==(B,A)
 */
static int compareName(char *nameA1, char *nameA2, 
                const char *nameB1, const char *nameB2) {

  if ( ((strcmp(nameA1,nameB1)==0) && (strcmp(nameA2,nameB2)==0)) ||
       ((strcmp(nameA1,nameB2)==0) && (strcmp(nameA2,nameB1)==0)) )
    return 0;
  else
    return 1;
}

/* compare(a,b)
 * Compare characters a and b, for use with qsort
 */
static int compareChar(const void *a, const void *b) {
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

/* atommap::SetDebug()
 * Set atommap debug
 */
void atommap::SetDebug(int debugIn) {
  debug=debugIn;
}

/* atommap::atomID()
 * Return the atomID of the given atom.
 */
const char *atommap::atomID(int atom) {
  if (atom<0 || atom>=natom) return NULL;
  return (M[atom].atomID);
}

/* atommap::Aname()
 * Return the parm atom name of the given atom.
 */
const char *atommap::Aname(int atom) {
  if (atom<0 || atom>=natom) return NULL;
  return (P->names[atom]);
}

/* atommap::getCut() 
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
    else if (strcmp(atom1,"S")==0) cut=2.05; // S-S gas-phase value; S=S is 1.49
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

/* atommap::calcDist()
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

  // Search for chiral centers by number of bonds
  for (i=0; i<natom; i++) {
    if (M[i].nbond==4) {
      // If >=3 bonds to single atoms, not chiral (e.g. -CH3)
      nbondj=0; // Count # bonds to single atoms
      for (j=0; j<M[i].nbond; j++) {
        nbondi = M[i].bond[j]; // nbondi is bonded to atom i
        if (M[nbondi].nbond==1) nbondj++;
      }
      if (nbondj<3) M[i].isChiral=true;
    }
  }

  // DEBUG - print bonding information
  if (debug>0) {
    mprintf("atommap: Atom Bond information.\n");
    for (i=0; i<natom; i++) {
      mprintf("  Atom %s_%i has %i bonds.",names[i],i,M[i].nbond);
      if (M[i].isChiral) mprintf(" CHIRAL!");
      mprintf("\n");
      for (j=0; j<M[i].nbond; j++) {
        nbondj=M[i].bond[j];
        mprintf("    to %s_%i\n",names[nbondj],nbondj);
      }
    }
  }

  return 0;
}

/* atommap::markAtomComplete()
 * If atom is mapped and all bonded atoms are mapped mark atom as completely 
 * mapped.
 * If printAtoms is true print isMapped value for this atom and all atoms
 * bonded to it.
 */
void atommap::markAtomComplete(int atom, bool printAtoms) {
  int nunique,bondatom;

  if (atom<0 || atom>=natom) return;
  if (!M[atom].isMapped && !printAtoms) return;
  if ( M[atom].complete && !printAtoms) return;
  nunique=0;
  for (int bond=0; bond < M[atom].nbond; bond++) {
    bondatom = M[atom].bond[bond];
    if (M[bondatom].isMapped) nunique++;
  }
  if (M[atom].isMapped && nunique==M[atom].nbond) {
    M[atom].complete=true;
  }
  if (printAtoms) {
    mprintf("  Atom %4i: %s-%1i |",atom,names[atom],(int)M[atom].isMapped);
    for (int bond=0; bond < M[atom].nbond; bond++) {
      bondatom = M[atom].bond[bond];
      mprintf(" %4i:%s-%1i",bondatom,names[bondatom],(int)M[bondatom].isMapped);
    }
    if (M[atom].complete)
      mprintf(" Atom is completely mapped.");
    mprintf("\n");
  }
}

/* atommap::markComplete()
 * Go through each atom in the map. If the atom is unique and all bonded
 * atoms are unique mark the atom as completely mapped.
 * Print bond information for each atom in the map, indicate whether
 * atoms are unique or not.
 */
void atommap::markComplete() {
  bool printAtoms = (debug>0);
  for (int atom=0; atom<natom; atom++) 
    markAtomComplete(atom, printAtoms);
}

/* atommap::determineAtomID()
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
        M[i].isUnique=false;
        M[j].isUnique=false;
        M[i].Nduplicated++;
        M[j].Nduplicated++;
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
      mprintf("  Atom %6i [%3i]: %s",i,M[i].Nduplicated,M[i].unique);
      if (M[i].isUnique) mprintf(" UNIQUE!");
      mprintf("\n");
    }
  }
}

/* atommap::BondIsRepeated()
 * Check if the atomID of the atom in bond <bond> bonded to atom <atom>
 * is the same as the atomID of any other non-mapped atom bonded to 
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
    // If bondedAtom2 is already mapped dont check it
    if (M[bondedAtom2].isMapped) continue;
    if (strcmp(atomID(bondedAtom),atomID(bondedAtom2))==0) return true;
  }
  return false;
}

/* atommap::setup()
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
    M[atom].isChiral=false;
    memset(M[atom].atomID,' ',ATOMIDLENGTH);
    memset(M[atom].unique,' ',UNIQUELENGTH);
    M[atom].isUnique=true; // Assume unique until proven otherwise
    M[atom].Nduplicated=0;
    M[atom].isMapped=false;
  }
  return 0;
}

// DEBUG
/* atommap::WriteMol2()
 * Write atommap out as a mol2 file, useful for checking bond info
 */
void atommap::WriteMol2(char *m2filename) {
  TrajectoryFile outfile;
  AtomMask M1;
  // Temporary parm to play with
  AmberParm *tmpParm;
  Frame *tmpFrame;
  ArgList tmpArg;

  // Create mask containing all atoms
  //for (int atom=0; atom<natom; atom++) Selected[atom]=atom;
  // Fake strip, just use as crap way to copy
  //tmpParm = P->modifyStateByMask(Selected,natom);
  //free(Selected);
  // Modify the bonds array to include this info
  //tmpParm->ResetBondInfo();
  //for (int atom=0; atom<natom; atom++) 
  //  for (int bond=0; bond < M[atom].nbond; bond++) 
  //    tmpParm->AddBond(atom, M[atom].bond[bond], 0);
  // Create mask with all mapped atoms
  for (int atom=0; atom<natom; atom++) {if (M[atom].isMapped) M1.AddAtom(atom);}
  // Strip so only mapped atoms remain
  tmpParm = P->modifyStateByMask(M1.Selected,M1.Nselected);
  tmpFrame = new Frame(M1.Nselected,NULL);
  tmpFrame->SetFrameFromMask(F, &M1);

  // Trajectory Setup
  tmpArg.Add((char*)"title\0");
  tmpArg.Add(m2filename);
  tmpArg.ResetAll();
  outfile.SetDebug(debug);
  outfile.SetupWrite(m2filename,&tmpArg,tmpParm,MOL2FILE);
  outfile.WriteFrame(0,tmpParm,tmpFrame->X,NULL,tmpFrame->box,tmpFrame->T);
  outfile.EndTraj();

  delete tmpParm;
  delete tmpFrame;
}
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

/* AtomMap::mapBondsToUnique()
 * For each atom R in reference already mapped to unique atom T in 
 * target, try to match up non-mapped reference atoms r bonded to R to 
 * non-mapped target atoms t bonded to T. Repeat until no more reference 
 * atoms can be mapped in this way.
 * Checking is very strict in this routine; r and t must be the only 
 * possible match and the atomIDs must match.
 * Return the total number of atoms mapped.
 */
int AtomMap::mapBondsToUnique(atommap *Ref, atommap *Tgt) {
  int atom,bond,r;
  int tatom,tbond,t;
  int numMappedAtoms=0;
  bool newSingle=true;

  while (newSingle) {

    // newSingle will be set back to true if any atoms are mapped
    newSingle=false;

    for (atom=0; atom < Ref->natom; atom++) {
      // Skip non-mapped atoms in Ref
      if (!Ref->M[atom].isMapped) continue;
      tatom = AMap[atom];
      // Check if map value is valid
      if (tatom<0) {
        mprintf("      Error: mapBondsToUnique: Ref %i:%s map value is invalid.\n",
                atom,Ref->Aname(atom));
        return -1;
      }
      // For each non-mapped atom bonded to Reference atom, try to 
      // find a matching non-mapped atom in Target by virtue of it being the
      // only possible match.
      for (bond=0; bond < Ref->M[atom].nbond; bond++) {
        r = Ref->M[atom].bond[bond];
        // Check that bonded atom r is not already mapped
        if (Ref->M[r].isMapped) continue;
        //mprintf("        Ref: Checking non-mapped %i:%s bonded to %i:%s\n",r,Ref->names[r],
        //        atom,Ref->names[atom]);
        // Check that non-mapped bonded ref atom r atomID is not the same as any 
        // other non-mapped bonded atomID.
        if ( Ref->BondIsRepeated(atom, bond) ) continue;
        // At this point r is the only one of its kind bonded to atom.
        // Check if there is an analogous atom bonded to unique Target atom
        // tatom.
        for (tbond=0; tbond < Tgt->M[tatom].nbond; tbond++) {
          t = Tgt->M[tatom].bond[tbond];
          // Check that bonded atom t is not already mapped 
          if (Tgt->M[t].isMapped) continue;
          //mprintf("          Tgt: Checking non-mapped %i:%s bonded to %i:%s\n",t,Tgt->names[t],
          //        tatom,Tgt->names[tatom]);
          // Check that non-mapped bonded tgt atom t atomID is not the same as 
          // any other non-mapped bonded atomID.
          if ( Tgt->BondIsRepeated(tatom, tbond) ) continue;
          // At this point t is the only one of its kind bonded to tatom.
          // Check if its atomID matches r. If so, map it.
          if ( strcmp(Ref->atomID(r), Tgt->atomID(t))==0 ) {
            if (debug>0) 
              mprintf("    Mapping tgt %i:%s to ref %i:%s based on single bond to unique.\n",
                      t,Tgt->Aname(t),r,Ref->Aname(r));
            AMap[r]=t;
            Ref->M[r].isMapped=true;
            Tgt->M[t].isMapped=true;
            newSingle=true;
            numMappedAtoms++;
          }
        } // End loop over atoms bonded to tatom
      } // End loop over atoms bonded to atom
      // Check if atom is completely mapped now
      Ref->markAtomComplete(atom,false);
      Tgt->markAtomComplete(tatom,false);
    } // End loop over ref atoms

  } // End loop over newSingle
  return numMappedAtoms;
}        

/* AtomMap::mapChiral()
 * Given two atommaps and a map relating the two, find chiral centers for
 * which at least 3 of the atoms have been mapped. Assign the remaining
 * two atoms based on improper dihedrals. 
 * Return the total number of mapped atoms.
 * NOTE: ONLY WORKS FOR SP3
 */
int AtomMap::mapChiral(atommap *Ref, atommap *Tgt) {
  int atom,tatom,bond,nunique,notunique_r,notunique_t;
  int uR[5], uT[5], r, t, nR[4], nT[4];
  double dR[4], dT[4], delta;
  int numMappedAtoms=0;

  for (atom=0; atom<Ref->natom; atom++) {
    // Skip non-mapped atoms
    if (!Ref->M[atom].isMapped) continue;
    //mprintf("DBG: mapChiral: Ref atom %i:%s\n",atom,Ref->P->names[atom]);
    tatom = AMap[atom];
    // Check that map value is valid
    if (tatom<0) {
      mprintf("      Error: mapChiral: Ref atom %i:%s map value is invalid.\n",
              atom,Ref->Aname(atom));
      return -1;
    }
    // If this Ref atom already completely mapped, skip
    if (Ref->M[atom].complete) {
      // Sanity check - if Ref atom is completely mapped, target should be
      // unless # atoms in Tgt and Ref are different.
      if (!Tgt->M[tatom].complete) {
        mprintf("Warning: mapChiral: Ref atom %i:%s is complete but Tgt atom %i:%s is not.\n",
                atom,Ref->Aname(atom),tatom,Tgt->Aname(tatom));
        //return 1;
      }
      continue;
    }
    // Check if this is a chiral center
    if (!Ref->M[atom].isChiral) continue;
    // If target atom is not a chiral center (e.g. due to diff # atoms)
    // mapping by chirality is not important for this reference, let
    // mapByIndex handle it.
    if (!Tgt->M[tatom].isChiral) {
      mprintf("Warning: mapChiral: Ref atom %i:%s is chiral but Tgt atom %i:%s is not!\n",
              atom,Ref->Aname(atom),tatom,Tgt->Aname(tatom));
      mprintf("         Marking Ref atom as non-chiral to try and map Tgt.\n");
      Ref->M[atom].isChiral=false;
      continue;
    }
    // Both atoms are chiral centers. Place bonded atoms (starting with 
    // central atom) in R and T.
    uR[0] = atom;
    uT[0] = tatom;
    nunique=1;
    notunique_r=0;
    // Look for mapped bonded ref and target atoms, and nonmapped reference atoms
    for (bond=0; bond<Ref->M[atom].nbond; bond++) {
      r = Ref->M[atom].bond[bond];
      t = AMap[r];
      // Bonded atom r is not mapped 
      if (!Ref->M[r].isMapped) {
        nR[notunique_r++] = r;
      // Bonded atom r is mapped. If a target was mapped to it
      // (i.e. it is the same atom) store it.
      } else {
        if (t>=0) {
          if (Ref->M[r].isMapped && Tgt->M[t].isMapped) {
            uR[nunique] = r;
            uT[nunique] = t;
            nunique++;
          }
        } 
      }
    }
    // Fill nT with nonmapped atoms from target
    notunique_t=0;
    for (bond=0; bond<Tgt->M[tatom].nbond; bond++) {
      t = Tgt->M[tatom].bond[bond];
      if (!Tgt->M[t].isMapped) nT[notunique_t++] = t;
    }
    // notunique_r may not be the same as notunique_t if the # atoms is different
    if (notunique_r!=notunique_t) 
      mprintf("Warning: Ref and Tgt do not have the same # of nonmapped atoms.\n");
    if (debug>0) { 
      mprintf("  Potential Chiral center %i_%s/%i_%s: Mapped atoms=%i, non-Mapped=%i/%i\n",
              atom,Ref->names[atom],tatom,Tgt->names[tatom],
              nunique,notunique_r,notunique_t);
      for (r=0; r<nunique; r++)
        mprintf("\t   Mapped\t%4i %4i\n",uR[r],uT[r]);
      for (r=0; r<notunique_r; r++)
        mprintf("\tNotMappedRef\t%4i\n",nR[r]);
      for (r=0; r<notunique_t; r++)
        mprintf("\tNotMappedTgt\t     %4i\n",nT[r]);
    }
    // If all atoms are unique no need to map
    // NOTE: Should be handled by complete check above.
    //if (nunique==5) continue;
    // Require at least 3 unique atoms for dihedral calc. 
    if (nunique<3) {
      if (debug>0) 
        mprintf("    Warning: Center has < 3 mapped atoms, dihedral cannot be calcd.\n");
      continue;
    }
    // Calculate reference improper dihedrals
    for (r=0; r<notunique_r; r++) {
      dR[r] = Torsion(Ref->F->Coord(uR[0]),Ref->F->Coord(uR[1]),
                      Ref->F->Coord(uR[2]),Ref->F->Coord(nR[r]));
      if (debug>1) mprintf("    Ref Improper %i [%3i,%3i,%3i,%3i]= %lf\n",r,
                           uR[0],uR[1],uR[2],nR[r],dR[r]);
    }
    // Calculate target improper dihedrals
    for (r=0; r<notunique_t; r++) {
      dT[r] = Torsion(Tgt->F->Coord(uT[0]),Tgt->F->Coord(uT[1]),
                      Tgt->F->Coord(uT[2]),Tgt->F->Coord(nT[r]));
      if (debug>1) mprintf("    Tgt Improper %i [%3i,%3i,%3i,%3i]= %lf\n",r,
                           uR[0],uR[1],uR[2],nT[r],dT[r]);
    }
    // Match impropers to each other using a cutoff
    // NOTE: 10.0 seems reasonable? Also there is currently no check for 
    //       repeated deltas.
    for (r=0; r<notunique_r; r++) {
      for (t=0; t<notunique_t; t++) {
        delta = dR[r] - dT[t];
        if (delta<0.0) delta=-delta;
        if (delta<10.0) {
          if (debug>0)
            mprintf("    Mapping tgt atom %i:%s to ref atom %i:%s based on chirality.\n",
                    nT[t],Tgt->Aname(nT[t]),nR[r],Ref->Aname(nR[r]) );
          AMap[ nR[r] ]=nT[t];
          numMappedAtoms++;
          // Once an atom has been mapped set its mapped flag
          Ref->M[nR[r]].isMapped=true;
          Tgt->M[nT[t]].isMapped=true;
        }
      }
    }
    // Check if ref atom or tgt atom is now completely mapped
    Ref->markAtomComplete(atom,false);
    Tgt->markAtomComplete(tatom,false);
  } // End loop over natom

  return numMappedAtoms;
}

/* AtomMap::mapUniqueRefToTgt()
 * If the number of atoms in Ref is different from Tgt, it is possible that
 * Tgt is missing atoms (or maybe vice versa). If the difference is not too
 * great it may be possible to look for an unmapped atom in Ref that has
 * same name and at least 1 matching bond (# bonds may be diff due to # atoms
 * so atomID cannot be used).
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
      if (debug>1) mprintf("        Attempting match of Tgt %i:%s to Ref %i:%s\n",
              t,Tgt->Aname(t),atom,Ref->Aname(atom));
      // Check that at least 1 bond is in common
      commonBond=0;
      for (bond=0; bond < Ref->M[atom].nbond; bond++) {
        // Check Map for ref bonded atom
        r = AMap[ Ref->M[atom].bond[bond] ];
        // If no mapping exists cant check it
        if (r<0) continue;
        if (debug>1) 
          mprintf("          Ref %i:%s bonded to %i:%s (%i:%s in tgt)\n",
                  atom, Ref->Aname(atom), 
                  Ref->M[atom].bond[bond], Ref->Aname( Ref->M[atom].bond[bond] ), 
                  r,Tgt->Aname(r));
        for (tbond=0; tbond < Tgt->M[t].nbond; tbond++) {
          if (debug>1)
            mprintf("            Tgt %i:%s bonded to %i:%s\n",t,Tgt->Aname(t),
                    Tgt->M[t].bond[tbond], Tgt->Aname(Tgt->M[t].bond[tbond]) );
          if (r == Tgt->M[t].bond[tbond]) 
            commonBond++;
        }
      }
      if (commonBond==0) continue;
      // This Tgt Name matches and at least 1 bond in common with Ref atom
      // Check that a match has not yet been found for ref
      if (match!=-1) {
        mprintf("      Warning: mapUniqueRefToTgt: Ref %i:%s has multiple potential matches\n",
                atom,Ref->Aname(atom));
        mprintf("               among Tgt [%i:%s, %i:%s]\n",
                t,Tgt->Aname(t),match,Tgt->Aname(match));
        return 0;
      }
      match = t;
    }
  }
  if (match==-1) return 0;
  // Only one match found - map it
  if (debug>0) 
    mprintf("    Mapping target %i:%s to unique ref %i:%s\n",match,Tgt->Aname(match),
            atom,Ref->Aname(atom));
  AMap[atom]=match;
  Ref->M[atom].isMapped=true;
  Tgt->M[match].isMapped=true;
  return 1;
}

/* AtomMap::mapByIndex()
 * Given to atommaps and a map relating the two, attempt to map any remaining
 * incomplete atoms by assuming the atom indices in reference and target are
 * in similar orders. At this point all unique atoms should have been mapped.
 * First, for each reference atom R check if R is unique but not mapped and 
 * attempt to match it to a non-mapped target based on local bonding 
 * environment (mapUniqueRefToTgt). Lastly, for reference atom R mapped to 
 * target atom T, compare the non-mapped atoms bonded to R (r) to the 
 * non-mapped atoms bonded to T (t). If the unique IDs of r and t match, map 
 * them. Otherwise if there is only one potential match between r and t map 
 * them.
 * Return the number of atoms mapped this way. 
 */
int AtomMap::mapByIndex(atommap *Ref, atommap *Tgt) {
  int atom,tatom,bond,tbond,numAtomsMapped,r,t;
  int match;

  numAtomsMapped=0;
  for (atom=0; atom<Ref->natom; atom++) {
    tatom = AMap[atom];
    // Check if no mapping exists for this atom 
    if (tatom<0) {
      // Check if reference atom is unique, but hasnt had a target mapped to it.
      // This can arise when the number of atoms in ref and tgt not equal.
      // If the difference in atoms is not too great (probably ~1) attempt
      // to look for a similar atom in tgt (name and index) and map it.
      if (Ref->M[atom].isUnique) {
        mprintf("      Warning: mapByIndex: Atom %i:%s in reference is unique but not mapped!\n",
                atom,Ref->Aname(atom));
        if (mapUniqueRefToTgt(Ref,Tgt,atom)) numAtomsMapped++;
      }
      continue;
    }
    // Skip over non-mapped atoms
    //if (!Ref->M[atom].isMapped) continue;

    // Check that num bonds match in Ref and target.
    // The # of bonds might not be equal if the # atoms in ref and tgt
    // not equal.
    if (Ref->M[atom].nbond!=Tgt->M[tatom].nbond) {
      mprintf("      Warning: mapByIndex: Ref atom %i:%s #bonds (%i) does not match Tgt atom %i:%s (%i)\n",
              atom,Ref->Aname(atom),Ref->M[atom].nbond,tatom,Tgt->Aname(tatom),Tgt->M[tatom].nbond);
      //return 1;
    }
    // Skip completely mapped atoms - check that both Ref and Tgt are complete
    // NOTE: This is ok if #atoms in ref > #atoms in tgt but not the other way around.
    if (Ref->M[atom].complete) {
      if (!Tgt->M[tatom].complete) {
        mprintf("Error: mapByIndex: Ref atom %i:%s is complete but Tgt atom %i:%s is not.\n",
                atom,Ref->Aname(atom),tatom,Tgt->Aname(tatom) );
        //return 1;
      }
      continue;
    }
    // This atom is mapped, but bonded atoms are not completely mapped. Try
    // to map the unmapped reference atoms bonded to <atom> to the unmapped
    // target atoms bonded to <tatom>. 
    if (debug>1)
      mprintf("DBG: Checking bonds of mapped Ref %i:%s against mapped Tgt %i:%s\n",
              atom,Ref->Aname(atom),tatom,Tgt->Aname(tatom));
    for (bond=0; bond < Ref->M[atom].nbond; bond++) {
      r = Ref->M[atom].bond[bond];
      if (debug>1) mprintf("\t\tRefBond %i:%s [%1i]\n",r,Ref->Aname(r),(int)Ref->M[r].isMapped);
      if (Ref->M[r].isMapped) continue;
      // Dont map atoms that are single-bonded to chiral centers; let
      // mapChiral take care of them.
      if (Ref->M[atom].isChiral && Ref->M[r].nbond==1) continue;
      match = -1;
      for (tbond=0; tbond < Tgt->M[tatom].nbond; tbond++) {
        t = Tgt->M[tatom].bond[tbond];
        if (debug>1) mprintf("\t\t\tTgtBond %i:%s [%1i]\n",t,Tgt->Aname(t),(int)Tgt->M[t].isMapped);
        if (Tgt->M[t].isMapped) continue;
        // Atom r bonded to atom, atom t bonded to tatom. r and t are not
        // yet mapped. Check if names match
        if (strcmp(Ref->names[r],Tgt->names[t])!=0) continue;
        // If the uniqueIDs of bonded atom r and bonded atom t match, map them now
        // NOTE: Scan for repeats?
        if (strcmp(Ref->M[r].unique,Tgt->M[t].unique)==0) {
          match = t;
          break;
        }
        // Store this atom t bonded to tatom as a potential match. If another
        // match has already been stored we cant tell these apart yet so ignore.
        if (match==-1) {
          match = t;
        } else {
          mprintf("      Warning: mapByIndex: Atom %i:%s bonded to Ref %i:%s has too many matches.\n",
                  r,Ref->Aname(r),atom,Ref->Aname(atom));
          match = -1;
          break;
        }
      } // End loop tbond over bonds in target atom
      // If a match was found, Map it
      if (match!=-1) {
        if (debug>0) mprintf("    Mapping Tgt %i:%s to Ref %i:%s based on name/bonding.\n",
                             match,Tgt->Aname(match),r,Ref->Aname(r));
        AMap[r] = match;
        Ref->M[r].isMapped=true;
        Tgt->M[match].isMapped=true;
        numAtomsMapped++;
      }
    } // End loop over atoms bonded to Ref atom
    // Check if atom is completely mapped now
    Ref->markAtomComplete(atom,false);
    Tgt->markAtomComplete(tatom,false);
  } // End loop over atoms

  return numAtomsMapped;
}


/* AtomMap::MapUniqueAtoms()
 * Map unique atoms in reference to unique atoms in target. If no atoms
 * can be mapped in this way, attempt to guess a starting point based
 * first on uniqueID, then by chirality.
 * Return number of atoms mapped.
 */
int AtomMap::MapUniqueAtoms(atommap *Ref, atommap *Tgt) {
  int refatom,targetatom;
  std::list<int> refGuess;
  std::list<int> tgtGuess;
  int numAtomsMapped=0;
  // Atoms have now been assigned IDs. Match up the unique strings in Ref with 
  // unique strings in target.
  for (refatom=0; refatom<Ref->natom; refatom++) {
    AMap[refatom]=-1;
    // If the ID of this reference atom is unique, look for same ID in target
    if (Ref->M[refatom].isUnique) {
      for (targetatom=0; targetatom<Tgt->natom; targetatom++) {
        // If ID of thie target atom is unique, check if it matches reference atom ID
        if (Tgt->M[targetatom].isUnique) {
          if ( strcmp(Tgt->M[targetatom].unique, Ref->M[refatom].unique)==0 ) {
            // Check that number of bonds is consistent
            if (Ref->M[refatom].nbond!=Tgt->M[targetatom].nbond) {
              mprintf("      Warning: AtomMap: Atoms R%i and T%i have same ID but different # bonds!\n",
                      refatom,targetatom);
            }
            AMap[refatom]=targetatom;
            Ref->M[refatom].isMapped=true;
            Tgt->M[targetatom].isMapped=true;
            numAtomsMapped++;
            if (debug>0) {
              mprintf("    Mapping Tgt %i:%s to Ref %i:%s based on unique ID\n",
                      targetatom,Tgt->Aname(targetatom),
                      refatom,Ref->Aname(refatom));
            }
          } // If unique strings match
        } // If target atom is unique
      } // Loop over target atoms
    } // If reference atom is unique
  } // Loop over reference atoms

  // If no unique atoms could be mapped it means the molecule is probably
  // very symmetric. At this point just try to guess a good starting
  // point. Map the first atoms that have a uniqueID duplicated only 1
  // time, preferably a chiral center.
  if (numAtomsMapped==0) {
    mprintf("      Warning: No unique atoms found, usually indicates highly symmetric system.\n");
    mprintf("               Trying to guess starting point.\n");
    for (refatom=0; refatom < Ref->natom; refatom++) {
      if (Ref->M[refatom].Nduplicated==1) {
        if (Ref->M[refatom].isChiral) 
          refGuess.push_front(refatom);
        else 
          refGuess.push_back(refatom);
      }
    }
    for (targetatom=0; targetatom < Tgt->natom; targetatom++) {
      if (Tgt->M[targetatom].Nduplicated==1) {
        if (Tgt->M[targetatom].isChiral)
          tgtGuess.push_front(targetatom);
        else
          tgtGuess.push_back(targetatom);
      }
    }
    if (refGuess.empty()) {
      mprintf("Error: AtomMap: Could not find starting point in reference.\n");
      return 0;
    }
    if (tgtGuess.empty()) {
      mprintf("Error: AtomMap: Could not find starting point in target.\n");
      return 0;
    }
    for (std::list<int>::iterator r=refGuess.begin(); r!=refGuess.end(); r++) {
      //mprintf("  Ref %i to ",*r);
      for (std::list<int>::iterator t=tgtGuess.begin(); t!=tgtGuess.end(); t++) {
        //mprintf("Tgt %i:",*t);
        if (strcmp(Ref->M[*r].unique, Tgt->M[*t].unique)==0) {
          //mprintf(" MATCH!\n");
          AMap[*r] = (*t);
          Ref->M[*r].isMapped=true;
          Tgt->M[*t].isMapped=true;
          numAtomsMapped++;
          mprintf("    Mapping Tgt %i:%s to Ref %i:%s based on guess.\n",*t,Tgt->Aname(*t),
                  *r,Ref->Aname(*r));
          break;
        }
        //mprintf("\n");
      }
      if (numAtomsMapped>0) break;
    }
    if (numAtomsMapped==0)
      mprintf("               Could not guess starting point.\n");
  } // End if numAtomsMapped==0
  return numAtomsMapped;
}

/* AtomMap::MapAtoms()
 * Map atoms in tgt to atoms in reference. First map any uniquely identified
 * atoms. Then map unmapped atoms that are the only one of their kind bonded 
 * to a unique or already mapped atom (mapBondsToUnique). Then map atoms based 
 * on chirality; if any atoms are mapped in this way check to see if 
 * mapBondsToUnique finds new atoms. Last try to guess mapping based on bonds 
 * (mapByIndex), which will also attempt to map atoms in Ref that are unique 
 * but not mapped to atoms in Tgt (which can happen e.g. if Tgt is missing 
 * atoms).
 * Negative return values from map... routines indicates error.
 */
int AtomMap::MapAtoms(atommap *Ref, atommap *Tgt) {
  bool mapatoms=true;
  int numAtomsMapped;
  int iterations=0;

  numAtomsMapped=MapUniqueAtoms(Ref, Tgt);
  // DEBUG
  //char name[1024];
  //sprintf(name,"Ref.%i.mol2",iterations);
  //Ref->WriteMol2(name);
  //sprintf(name,"Tgt.%i.mol2",iterations);
  //Tgt->WriteMol2(name);
  // END DEBUG
  if (debug>0)
    mprintf("*         MapUniqueAtoms: %i atoms mapped.\n",numAtomsMapped);
  if (numAtomsMapped==0) return 1;
  // Search for completely mapped atoms. If an atom and all atoms
  // it is bonded to are unique, mark the atom as completely mapped.
  RefMap.markComplete();
  TargetMap.markComplete();

  // Map remaining non-unique atoms
  while (mapatoms) {
    iterations++;
    // First assign based on bonds to unique (already mapped) atoms.
    numAtomsMapped=mapBondsToUnique(Ref,Tgt);
    // DEBUG
    //sprintf(name,"Ref.%i.u.mol2",iterations);
    //Ref->WriteMol2(name);
    //sprintf(name,"Tgt.%i.u.mol2",iterations);
    //Tgt->WriteMol2(name);
    // END DEBUG
    if (debug>0)
      mprintf("* [%3i] mapBondsToUnique: %i atoms mapped.\n",iterations,numAtomsMapped);
    if (numAtomsMapped<0) return 1;
    // Next assign based on chirality
    numAtomsMapped=mapChiral(Ref,Tgt);
    // DEBUG
    //sprintf(name,"Ref.%i.c.mol2",iterations);
    //Ref->WriteMol2(name);
    //sprintf(name,"Tgt.%i.c.mol2",iterations);
    //Tgt->WriteMol2(name);
    // END DEBUG
    if (debug>0)
      mprintf("* [%3i]        mapChiral: %i atoms mapped.\n",iterations,numAtomsMapped);
    if (numAtomsMapped<0) return 1;
    if (numAtomsMapped>0) continue;
    // Last assign based on index/element
    numAtomsMapped=mapByIndex(Ref,Tgt);
    // DEBUG
    //sprintf(name,"Ref.%i.i.mol2",iterations);
    //Ref->WriteMol2(name);
    //sprintf(name,"Tgt.%i.i.mol2",iterations);
    //Tgt->WriteMol2(name);
    // END DEBUG
    if (debug>0)
      mprintf("* [%3i]       mapByIndex: %i atoms mapped.\n",iterations,numAtomsMapped);
    if (numAtomsMapped<0) return 1;
    if (numAtomsMapped==0) mapatoms=false;
  }
  if (debug>0) mprintf("* %i iterations.\n",iterations);
  return 0;
}

/* AtomMap::init()
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

  // Check if number of atoms in each map is equal
  if (RefMap.natom!=TargetMap.natom) {
    mprintf("      AtomMap::init: Warning: # atoms in reference (%i) not equal\n",RefMap.natom);
    mprintf("                     to # atoms in target (%i).\n",TargetMap.natom);
  }

  // Allocate memory for atom map
  //   AMap[reference]=target
  AMap=(int*) malloc( RefMap.natom * sizeof(int));

  // Map atoms
  if (MapAtoms(&RefMap,&TargetMap)) return 1;

  // Print atom map and count # mapped atoms
  outputfile.SetupFile(outputname,WRITE,DATAFILE,UNKNOWN_TYPE,debug);
  outputfile.OpenFile();
  outputfile.IO->Printf("%-6s %4s %6s %4s\n","#TgtAt","Tgt","RefAt","Ref");
  for (refatom=0; refatom<RefMap.natom; refatom++) {
    targetatom=AMap[refatom];
    if (targetatom < 0) 
      outputfile.IO->Printf("%6s %4s %6i %4s\n","---","---",refatom+1,RefMap.Aname(refatom));
    else
      outputfile.IO->Printf("%6i %4s %6i %4s\n",targetatom+1,TargetMap.Aname(targetatom),
                            refatom+1,RefMap.Aname(refatom));
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
  if (maponly) return 0;

  // Check if not all atoms could be mapped
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
      // Since AMap[ ref ] = tgt but ref is now missing any stripped atoms,
      // the indices of AMap must be shifted to match
      refIndex=0; // The new index
      for (refatom=0; refatom<RefMap.natom; refatom++) {
        targetatom = AMap[refatom];
        if (targetatom<0)
          continue;
        else
          AMap[refIndex++]=targetatom;
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

/* AtomMap::setup()
 * If the current parm does not match the target parm, deactivate. Otherwise
 * replace current parm with mapped parm.
 */
int AtomMap::setup() {
  if (maponly) {
    mprintf("    ATOMMAP: maponly was specified, not using atom map during traj read.\n");
    return 0;
  }
  if (P->pindex!=TargetMap.P->pindex ||
      P->natom !=TargetMap.P->natom) 
  {
    mprintf("    ATOMMAP: Map for parm %s -> %s (%i atom).\n",TargetMap.P->parmName,
            RefMap.P->parmName,TargetMap.P->natom);
    mprintf("             Current parm %s (%i atom).\n",P->parmName,P->natom);
    mprintf("             Not using map for this parm.\n");
    return 1;
  }
  mprintf("    ATOMMAP: Map for parm %s -> %s (%i atom).\n",TargetMap.P->parmName,
          RefMap.P->parmName,TargetMap.P->natom);

  P = newParm;
  
  return 0;
}

/* AtomMap::action()
 * Modify the current frame based on the atom map. 
 */
int AtomMap::action() {
  if (maponly) return 0;
  for (int atom=0; atom < F->natom; atom++) 
    newFrame->SetCoord(atom, F->Coord(AMap[atom]));
  F = newFrame;
  return 0;
}

