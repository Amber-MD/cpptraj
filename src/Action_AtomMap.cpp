// AtomMap
#include <algorithm> //sort
#include <cstring> //memcpy
#include "Action_AtomMap.h"
#include "CpptrajStdio.h"
#include "TorsionRoutines.h"
#include "Bonds.h"
// DEBUG
#include "TrajectoryFile.h"
//#include <cstdio>

// atommap CONSTRUCTOR
atommap::atommap() {
  M=NULL;
  natom=0;
  names=NULL;
  mapFrame=NULL;
  mapParm=NULL;
  debug=0;
}

// atommap DESTRUCTOR
atommap::~atommap() {
  if (M!=NULL) delete[] M;
  if (names!=NULL) delete[] names; 
}

std::string atommap::NO_ID = "NULL";

// atommap::SetDebug()
/** Set atommap debug */
void atommap::SetDebug(int debugIn) {
  debug=debugIn;
}

// atommap::atomID()
/** Return the atomID of the given atom.  */
const std::string &atommap::atomID(int atom) {
  if (atom<0 || atom>=natom) return NO_ID;
  return (M[atom].atomID);
}

// atommap::Aname()
/** Return the parm atom name of the given atom.  */
const char *atommap::Aname(int atom) {
  return (mapParm->AtomName(atom));
}

// atommap::calcDist()
/** Determine which atoms are bonded to each other in a given set of atoms
  * based on how close they are and their identity.
  * Use bond info in parm. If no bond info generate it from the ref coord 
  * frame.
  * \return 1 on error, 0 on success.
  */
int atommap::calcDist() {
  std::vector<int> bondArray;

  // First check if bond info already present
  if ( mapParm->BondArray( bondArray ) ) {
    // No bond info present; create from frame
    mapParm->GetBondsFromCoords( mapFrame->X );
    // Try to get the bond info again.
    if ( mapParm->BondArray( bondArray ) ) {
      mprinterr("Error: AtomMap: Could not get bond information for %s\n",mapParm->parmName);
      return 1;
    }
  }
  // If here the bondArray was set up successfully. Array has format:
  //   Bond1:Atom1, Bond1:Atom2, Bond2:Atom1, Bond2:Atom2, Bond3:Atom1, ...
  // Fill in bond info for the map.
  for (std::vector<int>::iterator bondatom = bondArray.begin();
                                  bondatom != bondArray.end();
                                  bondatom++)
  {
    int atom1 = *bondatom;
    // Advance to atom2 in bond
    ++bondatom;
    int atom2 = *bondatom;
    // Bond atoms
    M[atom1].bond.push_back( atom2 );
    M[atom1].nbond++;
    M[atom2].bond.push_back( atom1 );
    M[atom2].nbond++;
    if (debug > 1) mprintf("\t%i:%s -- %i:%s\n",atom1+1, mapParm->AtomName(atom1),
                           atom2+1, mapParm->AtomName(atom2));
  }

  // Search for chiral centers by number of bonds
  for (int i=0; i<natom; i++) {
    // Sort the bonded atoms array by atom #
    sort( M[i].bond.begin(), M[i].bond.end() );
    if (M[i].nbond==4) {
      // If >=3 bonds to single atoms, not chiral (e.g. -CH3)
      int N_single_atoms=0; // Count # bonds to single atoms
      for (int idx=0; idx < M[i].nbond; idx++) {
        int bondedatom = M[i].bond[idx];
        if (M[bondedatom].nbond==1) ++N_single_atoms;
      }
      if (N_single_atoms<3) M[i].isChiral=true;
    }
  }

  // DEBUG - print bonding information
  if (debug>0) {
    mprintf("atommap: Atom Bond information.\n");
    for (int i=0; i<natom; i++) {
      mprintf("  Atom %s(%c)_%i has %i bonds.",mapParm->AtomName(i),names[i],i+1,M[i].nbond);
      if (M[i].isChiral) mprintf(" CHIRAL!");
      mprintf("\n");
      for (int j=0; j<M[i].nbond; j++) {
        int bondedatom=M[i].bond[j];
        mprintf("    to %s(%c)_%i\n",mapParm->AtomName(bondedatom),names[bondedatom],bondedatom+1);
      }
    }
  }

  return 0;
}

// atommap::markAtomComplete()
/** If atom is mapped and all bonded atoms are mapped mark atom as completely 
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
    mprintf("  Atom %4i: %c-%1i |",atom,names[atom],(int)M[atom].isMapped);
    for (int bond=0; bond < M[atom].nbond; bond++) {
      bondatom = M[atom].bond[bond];
      mprintf(" %4i:%c-%1i",bondatom,names[bondatom],(int)M[bondatom].isMapped);
    }
    if (M[atom].complete)
      mprintf(" Atom is completely mapped.");
    mprintf("\n");
  }
}

// atommap::markComplete()
/** Go through each atom in the map. If the atom is unique and all bonded
  * atoms are unique mark the atom as completely mapped.
  * Print bond information for each atom in the map, indicate whether
  * atoms are unique or not.
  */
void atommap::markComplete() {
  bool printAtoms = (debug>0);
  for (int atom=0; atom<natom; atom++) 
    markAtomComplete(atom, printAtoms);
}

// atommap::determineAtomID()
/** Give each atom an identifier based on what atoms are bonded to it. The
  * first part is the atom itself, followed by an alphabetized list of 
  * bonded atoms. So C in O=C-H2 would be CHHO.
  * Then determine which identifier strings are unique. 
  */
void atommap::determineAtomID() {
  // DEBUG
  //bool isRepeated;
  //int k, atom2;

  // Determine self id
  if (debug>0) mprintf("ATOM IDs:\n");
  for (int i=0; i<natom; i++) {
    M[i].atomID.clear();
    // Append the 1 char name of each bonded atom
    for (int j=0; j<M[i].nbond; j++) {
      int bondedatom=M[i].bond[j];
      M[i].atomID += names[bondedatom];
    }
    // Sort the bonded atom 1 char names
    // NOTE: Is it necessary for the atomID to be sorted, or just the Unique string?
    sort( M[i].atomID.begin(), M[i].atomID.end() );
    // Place this atoms 1 char name at the beginning
    M[i].atomID = names[i] + M[i].atomID;
    if (debug>0) mprintf("  Atom %i : %s\n",i,M[i].atomID.c_str());
  }

  // Create a unique ID for each atom based on Atom IDs
  for (int i=0; i<natom; i++) {
    M[i].unique = M[i].atomID;
    for (int j=0; j<M[i].nbond; j++) {
      int atom=M[i].bond[j];
      M[i].unique += M[atom].atomID;
    }
    sort( M[i].unique.begin(), M[i].unique.end() );
  }

  // Determine which unique IDs are duplicated - set isUnique flag
  for (int i=0; i<natom-1; i++) {
    for (int j=i+1; j<natom; j++) {
      if ( M[i].unique == M[j].unique ) {
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
    for (int i=0; i<natom; i++) {
      mprintf("  Atom %6i [%3i]: %s",i,M[i].Nduplicated,M[i].unique.c_str());
      if (M[i].isUnique) mprintf(" UNIQUE!");
      mprintf("\n");
    }
  }
}

// atommap::BondIsRepeated()
/** Check if the atomID of the atom in bond <bond> bonded to atom <atom>
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
    if ( atomID(bondedAtom) == atomID(bondedAtom2) ) return true;
  }
  return false;
}

// atommap::setup()
/** Allocate memory for atom map. In order to easily create the uniqueID 
  * strings the atom names need to be 1 char long. Convert chlorine 
  * to X for now, bromine to Y etc.
  */
int atommap::setup() {
  natom=mapParm->natom;
  names = new char[ natom ];
  // Set up 1 char atom names.
  for (int atom=0; atom<natom; atom++) {
    names[atom]=ConvertNameToChar(mapParm->AtomName(atom));
    // DEBUG
    if (debug>0) mprintf("  Atom %i:%s 1 char atom name: [%c]\n",
                         atom+1,mapParm->AtomName(atom),names[atom]);
  }
  // Allocate memory for atoms and initialize each atom
  M = new mapatom[ natom ];
  for (int atom=0; atom<natom; atom++) {
    //for (int bond=0; bond<MAXBONDS; bond++) M[atom].bond[bond]=-1;
    M[atom].nbond=0;
    M[atom].isChiral=false;
    M[atom].atomID.clear();
    M[atom].unique.clear();
    M[atom].isUnique=true; // Assume unique until proven otherwise
    M[atom].Nduplicated=0;
    M[atom].isMapped=false;
    M[atom].complete=false;
  }
  return 0;
}

// atommap::ResetMapping()
/** Reset all prior mapping-related information. */
void atommap::ResetMapping() {
  for (int atom=0; atom<natom; atom++) {
    M[atom].isMapped=false;
    M[atom].complete=false;
  }
}

// DEBUG
// atommap::WriteMol2()
/** Write atommap out as a mol2 file, useful for checking bond info
  */
void atommap::WriteMol2(char *m2filename) {
  TrajectoryFile outfile;
  AtomMask M1;
  // Temporary parm to play with
  AmberParm *tmpParm;
  Frame tmpFrame;
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
  tmpParm = mapParm->modifyStateByMask(M1.Selected,NULL);
  tmpFrame.SetupFrame(M1.Nselected,NULL);
  tmpFrame.SetFrameFromMask(mapFrame, &M1);

  // Trajectory Setup
  tmpArg.AddArg((char*)"title\0");
  tmpArg.AddArg(m2filename);
  outfile.SetDebug(debug);
  outfile.SetupWrite(m2filename,&tmpArg,tmpParm,MOL2FILE);
  outfile.WriteFrame(0,tmpParm,tmpFrame);
  outfile.EndTraj();

  delete tmpParm;
}
// ============================================================================

// CONSTRUCTOR
AtomMap::AtomMap() {
  AMap=NULL;
  newFrame=NULL;
  newParm=NULL;
  stripParm=NULL;
  maponly=false;
  rmsfit=false;
  rmsdata=NULL;
}

// DESTRUCTOR
AtomMap::~AtomMap() {
  if (AMap!=NULL) delete[] AMap;
  if (newFrame!=NULL) delete newFrame;
  if (newParm!=NULL) delete newParm;
  if (stripParm!=NULL) delete stripParm;
}

// AtomMap::mapBondsToUnique()
/** For each atom R in reference already mapped to unique atom T in 
  * target, try to match up non-mapped reference atoms r bonded to R to 
  * non-mapped target atoms t bonded to T. Repeat until no more reference 
  * atoms can be mapped in this way.
  * Checking is very strict in this routine; r and t must be the only 
  * possible match and the atomIDs must match.
  * \return the total number of atoms mapped.
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
          if ( Ref->atomID(r) == Tgt->atomID(t) ) {
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

// AtomMap::mapChiral()
/** Given two atommaps and a map relating the two, find chiral centers for
  * which at least 3 of the atoms have been mapped. Assign the remaining
  * two atoms based on improper dihedrals. 
  * \return the total number of mapped atoms.
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
      mprintf("  Potential Chiral center %i_%c/%i_%c: Mapped atoms=%i, non-Mapped=%i/%i\n",
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
      dR[r] = Torsion(Ref->mapFrame->Coord(uR[0]),Ref->mapFrame->Coord(uR[1]),
                      Ref->mapFrame->Coord(uR[2]),Ref->mapFrame->Coord(nR[r]));
      if (debug>1) mprintf("    Ref Improper %i [%3i,%3i,%3i,%3i]= %lf\n",r,
                           uR[0],uR[1],uR[2],nR[r],dR[r]);
    }
    // Calculate target improper dihedrals
    for (r=0; r<notunique_t; r++) {
      dT[r] = Torsion(Tgt->mapFrame->Coord(uT[0]),Tgt->mapFrame->Coord(uT[1]),
                      Tgt->mapFrame->Coord(uT[2]),Tgt->mapFrame->Coord(nT[r]));
      if (debug>1) mprintf("    Tgt Improper %i [%3i,%3i,%3i,%3i]= %lf\n",r,
                           uR[0],uR[1],uR[2],nT[r],dT[r]);
    }
    // Match impropers to each other using a cutoff. Note that all torsions
    // are in radians.
    // NOTE: 10.0 degrees seems reasonable? Also there is currently no 
    //       check for repeated deltas.
    for (r=0; r<notunique_r; r++) {
      for (t=0; t<notunique_t; t++) {
        delta = dR[r] - dT[t];
        if (delta<0.0) delta=-delta;
        if (delta<0.17453292519943295769236907684886) {
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

// AtomMap::mapUniqueRefToTgt()
/** If the number of atoms in Ref is different from Tgt, it is possible that
  * Tgt is missing atoms (or maybe vice versa). If the difference is not too
  * great it may be possible to look for an unmapped atom in Ref that has
  * same name and at least 1 matching bond (# bonds may be diff due to # atoms
  * so atomID cannot be used).
  * If only one name matches, probably safe to map it. 
  * \return 1 if the atom could be mapped, 0 otherwise.
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
    if ( Tgt->names[t] == Ref->names[atom] ) {
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

// AtomMap::mapByIndex()
/** Given to atommaps and a map relating the two, attempt to map any remaining
  * incomplete atoms by assuming the atom indices in reference and target are
  * in similar orders. At this point all unique atoms should have been mapped.
  * First, for each reference atom R check if R is unique but not mapped and 
  * attempt to match it to a non-mapped target based on local bonding 
  * environment (mapUniqueRefToTgt). Lastly, for reference atom R mapped to 
  * target atom T, compare the non-mapped atoms bonded to R (r) to the 
  * non-mapped atoms bonded to T (t). If the unique IDs of r and t match, map 
  * them. Otherwise if there is only one potential match between r and t map 
  * them.
  * \return the number of atoms mapped this way. 
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
      mprintf(
        "\tWarning: mapByIndex: Ref atom %i:%s #bonds (%i) does not match Tgt atom %i:%s (%i)\n",
        atom,Ref->Aname(atom),Ref->M[atom].nbond,tatom,Tgt->Aname(tatom),Tgt->M[tatom].nbond
      );
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
        if ( Ref->names[r] != Tgt->names[t] ) continue;
        // If the uniqueIDs of bonded atom r and bonded atom t match, map them now
        // NOTE: Scan for repeats?
        if ( Ref->M[r].unique == Tgt->M[t].unique ) {
          match = t;
          break;
        }
        // Store this atom t bonded to tatom as a potential match. If another
        // match has already been stored we cant tell these apart yet so ignore.
        if (match==-1) {
          match = t;
        } else {
          mprintf("\tWarning: mapByIndex: Atom %i:%s bonded to Ref %i:%s has too many matches.\n",
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

// AtomMap::MapUniqueAtoms()
/** Map unique atoms in reference to unique atoms in target. If no atoms
  * can be mapped in this way, attempt to guess a starting point based
  * first on uniqueID, then by chirality.
  * \return number of atoms mapped.
  */
int AtomMap::MapUniqueAtoms(atommap *Ref, atommap *Tgt) {
  int refatom,targetatom;
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
          if ( Tgt->M[targetatom].unique == Ref->M[refatom].unique ) {
            // Check that number of bonds is consistent
            if (Ref->M[refatom].nbond!=Tgt->M[targetatom].nbond) {
              mprintf("\tWarning: AtomMap: Atoms R%i and T%i have same ID but different # bonds!\n",
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

  return numAtomsMapped;
}

// AtomMap::MapWithNoUniqueAtoms()
/** If no unique atoms could be mapped it means the molecule is probably
  * very symmetric, so try to guess a good starting point. Map the first 
  * atoms that have a uniqueID duplicated only 1 time, preferably a chiral 
  * center. Try each pair of atoms and compute the resulting RMSD; use
  * the map with the lowest overall RMSD.
  * Note that the current implementation isn't very smart since it will 
  * try guess pairings that may have already been mapped in a previous
  * try.
  */
// NOTE: Also store the number of atoms mapped?
int AtomMap::MapWithNoUniqueAtoms( atommap *Ref, atommap *Tgt ) {
  std::list<int> refGuess;
  std::list<int> tgtGuess;
  double lowestRMS = 0;
  int *bestMap = NULL;
  int numAtomsMapped;
  double Rot[9], Trans[6];

  mprintf("      Warning: No unique atoms found, usually indicates highly symmetric system.\n");
  mprintf("               Trying to guess starting point.\n");
  //mprintf("DEBUG: Ref has %i atoms, Tgt has %i\n",Ref->natom, Tgt->natom);
  // Get a list of atoms in ref duplicated only once, preferably chiral
  for (int refatom=0; refatom < Ref->natom; refatom++) {
    if (Ref->M[refatom].Nduplicated==1) {
      if (Ref->M[refatom].isChiral) 
        refGuess.push_front(refatom);
      else 
        refGuess.push_back(refatom);
    }
  }
  if (refGuess.empty()) {
    mprintf("Error: AtomMap: Could not find starting point in reference.\n");
    return 1;
  }
  mprintf("Ref guess atoms:");
  for (std::list<int>::iterator r=refGuess.begin(); r!=refGuess.end(); r++)
    mprintf(" %i",(*r)+1);
  mprintf("\n");
  // Get a list of atoms in tgt duplicated only once, preferably chiral
  for (int targetatom=0; targetatom < Tgt->natom; targetatom++) {
    if (Tgt->M[targetatom].Nduplicated==1) {
      if (Tgt->M[targetatom].isChiral)
        tgtGuess.push_front(targetatom);
      else
        tgtGuess.push_back(targetatom);
    }
  }
  if (tgtGuess.empty()) {
    mprintf("Error: AtomMap: Could not find starting point in target.\n");
    return 1;
  }
  mprintf("Tgt guess atoms:");
  for (std::list<int>::iterator t=tgtGuess.begin(); t!=tgtGuess.end(); t++)
    mprintf(" %i",(*t)+1);
  mprintf("\n");
  // For each pair of atoms in refGuess and tgtGuess that have the same
  // ID string, guess that they are mapped and attempt to perform atom
  // mapping from there.
  for (std::list<int>::iterator r=refGuess.begin(); r!=refGuess.end(); r++) {
    for (std::list<int>::iterator t=tgtGuess.begin(); t!=tgtGuess.end(); t++) {
      if ( Ref->M[*r].unique == Tgt->M[*t].unique ) {
        // Reset any previous mapping
        for (int mapi=0; mapi < Ref->natom; mapi++) AMap[mapi]=-1;
        Ref->ResetMapping();
        Tgt->ResetMapping();
        //mprintf("  Ref %i (%i) to Tgt %i (%i) MATCH!\n",*r,Ref->natom,*t,Tgt->natom); // DEBUG
        // Map this guess
        AMap[(*r)] = (*t);
        Ref->M[(*r)].isMapped=true;
        Tgt->M[(*t)].isMapped=true;
        mprintf("    Mapping Tgt %i:%s to Ref %i:%s based on guess.\n",
                (*t)+1,Tgt->Aname(*t),
                (*r)+1,Ref->Aname(*r));
        // Attempt to complete mapping based on the guess
        if ( MapAtoms(Ref,Tgt) ) return 1;
        // Count number of mapped atoms
        numAtomsMapped=0;
        for (int mapi=0; mapi < Ref->natom; mapi++) if (AMap[mapi]!=-1) ++numAtomsMapped;
        // If < 3 atoms mapped this will cause a problem with RMSD
        if (numAtomsMapped<3) continue;
        // Score this mapping with an RMSD ---------------------------------
        // Set up a reference/target frame containing only mapped atoms
        int rmsIndex = 0;
        //mprintf("\tRMS fitting %i atoms from target to reference.\n",numAtomsMapped);
        rmsRefFrame.SetupFrame(numAtomsMapped,NULL);
        rmsTgtFrame.SetupFrame(numAtomsMapped,NULL);
        for (int refatom = 0; refatom < Ref->natom; refatom++) {
          int targetatom = AMap[refatom];
          if (targetatom!=-1) {
            rmsRefFrame.SetCoord(rmsIndex, Ref->mapFrame->Coord(refatom));
            rmsTgtFrame.SetCoord(rmsIndex, Tgt->mapFrame->Coord(targetatom));
            ++rmsIndex;
          }
        }
        double RmsVal = rmsTgtFrame.RMSD(&rmsRefFrame, Rot, Trans, false);
        mprintf("\tRMS fit (%i atoms) based on guess Tgt %i -> Ref %i, %lf\n",
                numAtomsMapped,(*t)+1, (*r)+1, RmsVal);
        // -----------------------------------------------------------------
        // If the current RmsVal is lower than the lowestRMS, store this map.
        if (bestMap==NULL || RmsVal < lowestRMS) {
          if (bestMap==NULL) bestMap = new int[ Ref->natom ];
          memcpy(bestMap, AMap, Ref->natom * sizeof(int));
          lowestRMS = RmsVal;
        }
      }
    } // End loop over tgt guesses
  } // End loop over ref guesses

  // If bestMap is NULL something went wrong. Otherwise set AMap to best map.
  if (bestMap==NULL) {
    mprinterr("Error: AtomMap::MapWithNoUniqueAtoms: Could not guess starting point.\n");
    return 1;
  } else {
    memcpy(AMap, bestMap, Ref->natom * sizeof(int));
    delete[] bestMap;
  }
  return 0;
}

// AtomMap::MapAtoms()
/** Map atoms in tgt to atoms in reference. Assumes that any uniquely 
  * identified atoms have already been mapped. First map unmapped atoms 
  * that are the only one of their kind bonded to a unique or already 
  * mapped atom (mapBondsToUnique). Then map atoms based on chirality; 
  * if any atoms are mapped in this way check to see if mapBondsToUnique 
  * finds new atoms. Last try to guess mapping based on bonds (mapByIndex), 
  * which will also attempt to map atoms in Ref that are unique but not 
  * mapped to atoms in Tgt (which can happen e.g. if Tgt is missing atoms).
  * Negative return values from mapXXX routines indicates error.
  * \return 0 on success, 1 on error.
  */
int AtomMap::MapAtoms(atommap *Ref, atommap *Tgt) {
  bool mapatoms=true;
  int numAtomsMapped;
  int iterations=0;

  // DEBUG
  //char name[1024];
  //sprintf(name,"Ref.%i.mol2",iterations);
  //Ref->WriteMol2(name);
  //sprintf(name,"Tgt.%i.mol2",iterations);
  //Tgt->WriteMol2(name);
  // END DEBUG
  // Search for completely mapped atoms. If an atom and all atoms
  // it is bonded to are unique, mark the atom as completely mapped.
  RefMap.markComplete();
  TargetMap.markComplete();

  // Map remaining non-unique atoms
  while (mapatoms) {
    ++iterations;
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

// AtomMap::init()
/** Expected call: atommap <target> <reference> [mapout <filename>] [maponly]
  *                        [rmsfit [ rmsout <rmsout> ]]
  * Attempt to create a map from atoms in target to atoms in reference solely
  * based on how they are bonded (not how they are named). 
  */
int AtomMap::init() {
  char *refName, *targetName, *outputname;
  char *rmsout = NULL;
  CpptrajFile outputfile;
  int refIndex, targetIndex;
  int refatom,targetatom;
  int numMappedAtoms;
  AtomMask *M1;
  
  RefMap.SetDebug(debug);
  TargetMap.SetDebug(debug);

  // Get Args
  outputname=actionArgs.getKeyString("mapout",NULL);
  maponly = actionArgs.hasKey("maponly");
  rmsfit = actionArgs.hasKey("rmsfit");
  if (rmsfit)
    rmsout = actionArgs.getKeyString("rmsout",NULL);

  targetName=actionArgs.getNextString();
  refName=actionArgs.getNextString();
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
  RefMap.mapFrame=FL->GetFrame(refIndex);
  // Get reference parm
  RefMap.mapParm=FL->GetFrameParm(refIndex);
  if (RefMap.mapFrame==NULL || RefMap.mapParm==NULL) {
    mprintf("AtomMap::init: Error: Could not get reference frame %s\n",refName);
    return 1;
  }
  // Get target index based on filename
  targetIndex=FL->GetFrameIndex(targetName);
  // Get target frame 
  TargetMap.mapFrame=FL->GetFrame(targetIndex);
  // Get target parm
  TargetMap.mapParm=FL->GetFrameParm(targetIndex);
  if (TargetMap.mapFrame==NULL || TargetMap.mapParm==NULL) {
    mprintf("AtomMap::init: Error: Could not get target frame %s\n",targetName);
    return 1;
  }

  mprintf("    ATOMMAP: Atoms in trajectories associated with parm %s will be\n",
          TargetMap.mapParm->parmName);
  mprintf("             mapped according to parm %s.\n",RefMap.mapParm->parmName);
  if (outputname!=NULL)
    mprintf("             Map will be written to %s\n",outputname);
  if (maponly)
    mprintf("             maponly: Map will only be written, not used in trajectory read.\n");
  if (!maponly && rmsfit) {
    mprintf("             rmsfit: Will rms fit mapped atoms in tgt to reference.\n");
    if (rmsout!=NULL) {
      rmsdata = DSL->Add(DOUBLE,actionArgs.getNextString(),"RMSD");
      if (rmsdata==NULL) return 1;
      DFL->Add(rmsout,rmsdata);
    }
  }

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
  AMap=new int[ RefMap.natom ]; 
  // Map unique atoms
  numMappedAtoms=MapUniqueAtoms(&RefMap, &TargetMap);
  if (debug>0)
    mprintf("*         MapUniqueAtoms: %i atoms mapped.\n",numMappedAtoms);
  // If no unique atoms mapped system is highly symmetric and needs to be
  // iteratively mapped. Otherwise just map remaining atoms.
  if (numMappedAtoms==0) { 
    if (MapWithNoUniqueAtoms(&RefMap,&TargetMap)) return 1;
  } else {
    if (MapAtoms(&RefMap,&TargetMap)) return 1;
  }

  // Print atom map and count # mapped atoms
  numMappedAtoms = 0;
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
      ++numMappedAtoms;
    } //else {
    //  mprintf("* Could not map any TargetAtom to RefAtom %6i(%4s)\n",
    //                  refatom,RefMap.P->names[refatom]);
    //}
  }
  outputfile.CloseFile();
  mprintf("      %i total atoms were mapped.\n",numMappedAtoms);
  if (maponly) return 0;

  // If rmsfit is specified, an rms fit of target to reference will be
  // performed using all atoms that were successfully mapped from 
  // target to reference.
  if (rmsfit) {
    int rmsRefIndex = 0;
    // Set up a reference frame containing only mapped reference atoms
    rmsRefFrame.SetupFrame(numMappedAtoms,NULL);
    for (refatom = 0; refatom < RefMap.natom; refatom++) {
      targetatom = AMap[refatom];
      if (targetatom!=-1) rmsRefFrame.SetCoord(rmsRefIndex++, RefMap.mapFrame->Coord(refatom));
    }
    // Prepare target frame to hold mapped atoms
    rmsTgtFrame.SetupFrame(numMappedAtoms,NULL);
    mprintf("      rmsfit: Will rms fit %i atoms from target to reference.\n",numMappedAtoms);
    return 0;
  }

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
      //stripParm = RefMap.mapParm->modifyStateByMask(M1->Selected, numMappedAtoms, NULL);
      stripParm = RefMap.mapParm->modifyStateByMask(M1->Selected, NULL);
      // Strip reference frame
      newFrame = new Frame();
      newFrame->SetupFrame(numMappedAtoms,RefMap.mapParm->mass);
      newFrame->SetFrameFromMask(RefMap.mapFrame, M1);
      delete M1;
      // Replace reference with stripped versions
      if (FL->ReplaceFrame(refIndex, newFrame, stripParm)) {
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
      mprintf("Warning: AtomMap: Not all atoms were mapped. Frames will not be modified.\n");
      maponly=true;
    }
  }

  if (!maponly) {
    // Set up new Frame
    newFrame = new Frame();
    newFrame->SetupFrame(TargetMap.natom,TargetMap.mapParm->mass);

    // Set up new Parm
    newParm = TargetMap.mapParm->modifyStateByMap(AMap);
  }

  return 0;
}

// AtomMap::setup()
/** If the current parm does not match the target parm, deactivate. Otherwise
  * replace current parm with mapped parm.
  */
int AtomMap::setup() {
  if (maponly) {
    mprintf("    ATOMMAP: maponly was specified, not using atom map during traj read.\n");
    return 0;
  }
  if (currentParm->pindex!=TargetMap.mapParm->pindex ||
      currentParm->natom !=TargetMap.mapParm->natom) 
  {
    mprintf("    ATOMMAP: Map for parm %s -> %s (%i atom).\n",TargetMap.mapParm->parmName,
            RefMap.mapParm->parmName,TargetMap.mapParm->natom);
    mprintf("             Current parm %s (%i atom).\n",currentParm->parmName,currentParm->natom);
    mprintf("             Not using map for this parm.\n");
    return 1;
  }
  if (rmsfit) {
    mprintf("    ATOMMAP: rmsfit specified, %i atoms.\n",rmsRefFrame.natom);
    return 0;
  }
  mprintf("    ATOMMAP: Map for parm %s -> %s (%i atom).\n",TargetMap.mapParm->parmName,
          RefMap.mapParm->parmName,TargetMap.mapParm->natom);

  currentParm = newParm;
  
  return 0;
}

// AtomMap::action()
/** Modify the current frame based on the atom map. 
  */
int AtomMap::action() {
  double Rot[9], Trans[6], R;
  if (maponly) return 0;

  // Perform RMS fit on mapped atoms only
  if (rmsfit) {
    int rmsTgtIndex = 0;
    for (int refatom = 0; refatom < RefMap.natom; refatom++) {
      int targetatom = AMap[refatom];
      if (targetatom!=-1) rmsTgtFrame.SetCoord(rmsTgtIndex++, currentFrame->Coord(targetatom));
    }
    R = rmsTgtFrame.RMSD(&rmsRefFrame, Rot, Trans, false);
    currentFrame->Trans_Rot_Trans(Trans,Rot);
    if (rmsdata!=NULL)
      rmsdata->Add(frameNum, &R);
    return 0;
  }

  // Modify the current frame
  for (int atom=0; atom < currentFrame->natom; atom++) 
    newFrame->SetCoord(atom, currentFrame->Coord(AMap[atom]));
  currentFrame = newFrame;
  return 0;
}

