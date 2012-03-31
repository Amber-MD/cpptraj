// AmberParm.cpp
#include <cstring>     // memcpy, strcpy
#include <cstdio>      // sprintf
#include "AmberParm.h" 
#include "PtrajMask.h" // parseMaskString
#include "Constants.h" // ELECTOAMBER
#include "CpptrajStdio.h"
// For searching for bonds by distance (PDB etc)
#include "DistRoutines.h"

// CONSTRUCTOR
AmberParm::AmberParm() {
  debug=0;
  parmfileName=NULL;
  parmName=NULL;
  pindex=0;
  parmFrames=0;
  parmCoords=NULL;

  NbondsWithH=0;
  NbondsWithoutH=0;
  bondsh=NULL;
  bonds=NULL;
  names=NULL;
  resnames=NULL;
  types=NULL;
  resnums=NULL;
  natom=0;
  nres=0;
  finalSoluteRes=0;
  molecules=0;
  firstSolvMol=-1;
  atomsPerMol=NULL;
  mass=NULL;
  charge=NULL;
  Box[0]=0.0; Box[1]=0.0; Box[2]=0.0;
  Box[3]=0.0; Box[4]=0.0; Box[5]=0.0;
  boxType=NOBOX;

  solventMask=NULL;
  solventMolecules=0;
  solventMoleculeStart=NULL;
  solventMoleculeStop=NULL;
  solventAtoms=0;

  SurfaceInfo=NULL;
  numSoluteAtoms=0;

  numex=NULL;
  atype_index=NULL;
  at_num=NULL;
  NB_index=NULL;
  LJ_A=NULL;
  LJ_B=NULL;
  excludedAtoms=NULL;
  radius_set=NULL;
  gb_radii=NULL;
  gb_screen=NULL;
  ntypes=0;
  nnb=0;

  bond_rk=NULL;
  bond_req=NULL;
  angle_tk=NULL;
  angle_teq=NULL;
  dihedral_pk=NULL;
  dihedral_pn=NULL;
  dihedral_phase=NULL;
  scee_scale=NULL;
  scnb_scale=NULL;
  solty=NULL;
  anglesh=NULL;
  angles=NULL;
  dihedralsh=NULL;
  dihedrals=NULL;
  asol=NULL;
  bsol=NULL;
  hbcut=NULL;
  itree=NULL;
  join_array=NULL;
  irotat=NULL;

  numbnd=0;
  numang=0;
  numdih=0;
  NanglesWithH=0;
  NanglesWithoutH=0;
  NdihedralsWithH=0;
  NdihedralsWithoutH=0;
  natyp=0;
  nphb=0;
}

// DESTRUCTOR
AmberParm::~AmberParm() {
  if (parmfileName!=NULL) delete[] parmfileName;
  if (parmName!=NULL) delete[] parmName;

  if (bondsh!=NULL) delete[] bondsh;
  if (bonds!=NULL) delete[] bonds;
  if (names!=NULL) delete[] names;
  if (resnames!=NULL) delete[] resnames;
  if (types!=NULL) delete[] types;
  if (resnums!=NULL) delete[] resnums;
  if (atomsPerMol!=NULL) delete[] atomsPerMol;
  if (mass!=NULL) delete[] mass;
  if (charge!=NULL) delete[] charge;

  if (solventMask!=NULL) delete[] solventMask;
  if (solventMoleculeStart!=NULL) delete[] solventMoleculeStart;
  if (solventMoleculeStop!=NULL) delete[] solventMoleculeStop;

  if (SurfaceInfo!=NULL) delete[] SurfaceInfo;
  if (parmCoords!=NULL) delete[] parmCoords;

  if (numex!=NULL) delete[] numex;
  if (atype_index!=NULL) delete[] atype_index;
  if (at_num!=NULL) delete[] at_num;
  if (NB_index!=NULL) delete[] NB_index;
  if (LJ_A!=NULL) delete[] LJ_A;
  if (LJ_B!=NULL) delete[] LJ_B;
  if (excludedAtoms!=NULL) delete[] excludedAtoms;
  if (radius_set!=NULL) delete[] radius_set;
  if (gb_radii!=NULL) delete[] gb_radii;
  if (gb_screen!=NULL) delete[] gb_screen;

  if (bond_rk!=NULL) delete[] bond_rk;
  if (bond_req!=NULL) delete[] bond_req;
  if (angle_tk!=NULL) delete[] angle_tk;
  if (angle_teq!=NULL) delete[] angle_teq;
  if (dihedral_pk!=NULL) delete[] dihedral_pk;
  if (dihedral_pn!=NULL) delete[] dihedral_pn;
  if (dihedral_phase!=NULL) delete[] dihedral_phase;
  if (scee_scale!=NULL) delete[] scee_scale;
  if (scnb_scale!=NULL) delete[] scnb_scale;
  if (solty!=NULL) delete[] solty;
  if (anglesh!=NULL) delete[] anglesh;
  if (angles!=NULL) delete[] angles;
  if (dihedralsh!=NULL) delete[] dihedralsh;
  if (dihedrals!=NULL) delete[] dihedrals;
  if (asol!=NULL) delete[] asol;
  if (bsol!=NULL) delete[] bsol;
  if (hbcut!=NULL) delete[] hbcut;
  if (itree!=NULL) delete[] itree;
  if (join_array!=NULL) delete[] join_array;
  if (irotat!=NULL) delete[] irotat;
}

// SetDebug()
/** Set the debug level.  */
void AmberParm::SetDebug(int debugIn) {
  debug = debugIn;
  if (debug>0) mprintf("AmberParm debug set to %i\n",debug);
}

// ---------- MASK PARSER ROUTINES --------------------------------------------- 
// AmberParm::SetupAtomMask()
/** Given an atommask which has a mask expression, set it up for this topology.
  * Most routines will not call this directly, but will call SetupIntegerMask
  * or SetupCharMask.
  * \param atommaskIn The AtomMask to set-up.
  * \param Xin Optional coordinates for setting up distance-based masks.
  * \param isCharMask True if atommaskIn should be set up as a character
  *                   mask, false if it should be set up as an integer mask.
  * \return 1 on error, 0 on success.
  */
int AmberParm::SetupAtomMask(AtomMask &atommaskIn, double *Xin, bool isCharMask) {
  char *atommask_postfix;
  char *ptraj_charmask;

  atommask_postfix = atommaskIn.PostfixExpression();
  if (atommask_postfix == NULL) {
    mprinterr("Error: AmberParm::SetupAtomMask: AtomMask Postfix expression is NULL.\n");
    return 1;
  }
  if (debug>1) mprintf("\tAmberParm::SetupAtomMask: Parsing postfix [%s]\n",atommask_postfix); 
 
  ptraj_charmask = parseMaskString(atommask_postfix, natom, nres, names, resnames,
                                   resnums, Xin, types, debug);

  if (ptraj_charmask==NULL) {
    mprinterr("    Error: AmberParm::SetupAtomMask: charmask is NULL.\n");
    return 1;
  }

  if (isCharMask)
    atommaskIn.SetupCharMask( ptraj_charmask, natom, debug);
  else
    atommaskIn.SetupMask( ptraj_charmask, natom, debug );

  // Free the character mask, no longer needed.
  delete[] ptraj_charmask;

  return 0;
} 

// AmberParm::SetupIntegerMask()
/** Set up the given AtomMask as an integer mask. */  
int AmberParm::SetupIntegerMask(AtomMask &atommaskIn, double *Xin) {
  return SetupAtomMask(atommaskIn, Xin, false);
}

// AmberParm::SetupCharMask()
/** Set up the given AtomMask as a character mask. */
int AmberParm::SetupCharMask(AtomMask &atommaskIn, double *Xin) {
  return SetupAtomMask(atommaskIn, Xin, true);
}

// AmberParm::NumMoleculesInMask()
/** Given an atom mask that has already been set up, return the number
  * of distinct molecules contained by that mask. Assumes the mask
  * is sorted.
  */
int AmberParm::NumMoleculesInMask(AtomMask &atommaskIn) {
  int numMolecules, lastMolecule;

  if (atomsPerMol == NULL) return 0;
  numMolecules = 0;
  lastMolecule = -1;
  for (std::vector<int>::iterator atom = atommaskIn.Selected.begin();
                                  atom != atommaskIn.Selected.end();
                                  atom++)
  {
    if (lastMolecule == -1) {
      lastMolecule = atomToMolecule( *atom );
      ++numMolecules;
    } else {
      int thisMolecule = atomToMolecule( *atom );
      if (thisMolecule != lastMolecule) {
        ++numMolecules;
        lastMolecule = thisMolecule;
      }
    }
  }
  return numMolecules;
}

// ---------- ROUTINES FOR ACCESSING RESIDUE AND ATOM INFORMATION -------------- 
// AmberParm::ResName()
/** Given a residue number, set buffer with residue name and number with format:
  * <resname[res]><res+1>, e.g. ARG_11. Replace any blanks in resname with '_'.
  */
void AmberParm::ResName(char *buffer, int res) {
  char rname[NAMESIZE];
  if (res<0 || res>=nres) return;
  rname[0]=resnames[res][0];
  rname[1]=resnames[res][1];
  rname[2]=resnames[res][2];
  if (resnames[res][3]==' ') 
    rname[3]='_';
  else
    rname[3]=resnames[res][3];
  rname[4]='\0';
  sprintf(buffer,"%s%i",rname,res+1);
}

// AmberParm::ResAtomName()
/** Given an atom number, set buffer with residue name and number along with
  * atom name with format: <resname[res]><res+1>@<atomname>, e.g. ARG_11@CA.
  * Replace any blanks in resname with '_'.
  */
void AmberParm::ResAtomName(char *buffer, int atom) {
  int res;
  char rname[NAMESIZE];
  if (atom<0 || atom>=natom) return;
  res = atomToResidue(atom);
  rname[0]=resnames[res][0];
  rname[1]=resnames[res][1];
  rname[2]=resnames[res][2];
  if (resnames[res][3]==' ') 
    rname[3]='_';
  else
    rname[3]=resnames[res][3];
  rname[4]='\0';
  sprintf(buffer,"%s%i@%s",rname,res+1,names[atom]);
}

// AmberParm::ResidueName()
/** Return pointer to name of given residue.  */
char *AmberParm::ResidueName(int res) {
  if (resnames==NULL) {
    mprintf("Internal Error: AmberParm::ResidueName: Residue names not set!\n");
    return NULL;
  }
  if (res>-1 && res<nres)
    return (char*)resnames[res];
  return NULL;
}

// AmberParm::FindAtomInResidue()
/** Find the atom # of the specified atom name in the given residue.
  * \param res Residue number to search.
  * \param atname Atom name to find.
  * \return the atom number of the specified atom if found in the given residue.
  * \return -1 if atom not found in given residue.
  */
int AmberParm::FindAtomInResidue(int res, char *atname) {
  if (res < 0 || res >= nres) return -1;
  for (int atnum = resnums[res]; atnum < resnums[res+1]; atnum++) {
    if (strcmp(names[atnum],atname)==0) return atnum;
  }
  return -1;
}

// AmberParm::ResAtomRange()
/** Set the first and last+1 atoms for the given residue */
int AmberParm::ResAtomRange(int resIn, int *startatom, int *stopatom) {
  if (resIn < 0 || resIn >= nres) {
    *startatom = 0;
    *stopatom = 0;
    return 1;
  }
  *startatom = resnums[resIn];
  *stopatom = resnums[resIn+1];
  return 0;
}

// AmberParm::FinalSoluteRes()
/** Return the last residue considered solute. */
int AmberParm::FinalSoluteRes() {
  if (hasSolventInfo)
    return finalSoluteRes;
  else
    return nres;
}

// AmberParm::AtomName()
/** Return pointer to name of given atom. */
char *AmberParm::AtomName(int atom) {
  if (atom<0 || atom >= natom) return NULL;
  return (char*)names[atom];
}

// AmberParm::AtomNameIs()
/** Return true if atom name matches input */
bool AmberParm::AtomNameIs(int atom, const char *nameIn) {
  if (atom<0 || atom >= natom) return false;
  if (strcmp(names[atom], nameIn)==0) return true;
  return false;
}

// AmberParm::AtomElementIs()
/** Return true if atom element matches input */
bool AmberParm::AtomElementIs(int atom, AtomicElementType elementIn) {
  if (atom<0 || atom >= natom) return false;
  AtomicElementType element = ElementFromName( names[atom] );
  if (element == elementIn) return true;
  return false;
}

// ---------- ROUTINES FOR ACCESSING INTERNAL PARAMETERS -----------------------
// AmberParm::SetupExcludedAtomsList()
/** Set up a vector that contains the exclusion list for each atom in the
  * given mask. The exclusion list for each atom will be terminated by a
  * -1 in order to easily check for the end of the excluded atoms list
  * which is useful e.g. when iterating over atom pairs in Action_Pairwise.
  * Atoms not in the mask will have a completely empty exclusion list.
  * There should be ((N^2 - N) / 2) - SUM(Numex) interactions, however
  * since the excluded atoms list can include 0 as a placeholder (thereby
  * having a Numex of 1) there may be more interactions than expected,
  * so subtract the actual number of excluded atoms from the total as
  * the overall exclusion list is built.
  * \param maskIn AtomMask (integer) containing atoms that will be calcd.
  * \param exclusionList Vector of vectors containing exclusion list for
  *                      each atom.
  * \return the total number of interactions between atoms in mask.
  * \return -1 on error.
  */
int AmberParm::SetupExcludedAtomsList(AtomMask &maskIn, 
                                      std::vector< std::vector<int> > &exclusionList) 
{
  int N_interactions;

  // Resize excluded atom list for current number of atoms
  exclusionList.clear();
  exclusionList.resize( natom );

  // Calculate the initial number of interactions
  N_interactions = (((maskIn.Nselected * maskIn.Nselected) - maskIn.Nselected) / 2);

  // If no exclusion information, just set -1 for each atom's exclusion list.
  // Since no atoms are excluded the total number of interactions do not need
  // to be modified.
  if (numex==NULL) {
    for (std::vector<int>::iterator maskatom1 = maskIn.Selected.begin(); 
                                    maskatom1 != maskIn.Selected.end(); 
                                    maskatom1++)
      exclusionList[ *maskatom1 ].push_back( -1 );

  // If there is exclusion information, for each atom in the mask add its 
  // excludedAtoms (that are also in the mask) and subtract the number excluded 
  // from the overall number of interactions.
  } else {
    // Need to ensure that atoms in the exclusion list are also in the mask, so set 
    // up a copy of maskIn that is a character mask for this purpose.
    AtomMask Mask2 = maskIn;
    if (Mask2.ConvertMaskType()) {
      mprinterr("Error: AmberParm::SetupExcludedAtomsList: Could not convert mask.\n");
      return -1;
    }
    int natex_idx = 0; // Index into excludedAtoms (NATEX)
    for (int atom = 0; atom < natom; atom++) {
      // Check if this atom is in the mask. If so, put all of its excluded
      // atoms that are also in the mask into the excluded list for this
      // atom. If not, just increment the index into excludedAtoms.
      if ( Mask2.AtomInCharMask( atom ) ) {
        for (int e_idx = 0; e_idx < numex[atom]; e_idx++) {
          int excluded_atom = excludedAtoms[ natex_idx++ ];
          if ( Mask2.AtomInCharMask( excluded_atom ) )
            exclusionList[ atom ].push_back( excluded_atom );
        }
        // Push back a final -1. Since no atom will ever have this number
        // it signals the end of the exclusion list for this atom. 
        exclusionList[ atom ].push_back( -1 );
        // The total number of excluded atoms for this atom is now the size
        // of the vector minus one (for the terminal -1).
        int total_excluded = (int)exclusionList[ atom ].size();
        N_interactions -= (total_excluded - 1);
      } else {
        natex_idx += numex[atom];
      }
    } // END loop over all atoms
  }
  // DEBUG - Print total number of interactions for this parm
  if (debug > 0) {
    mprintf("\t%i interactions for %s.\n",N_interactions,parmName);

    // DEBUG - Print exclusion list for each atom
    if (debug > 1) {
      for (unsigned int atom = 0; atom < exclusionList.size(); atom++) {
        mprintf("\t%8u:",atom + 1);
        for (std::vector<int>::iterator eat = exclusionList[atom].begin();
                                        eat != exclusionList[atom].end();
                                        eat++)
        {
          mprintf(" %i",*eat + 1);
        }
        mprintf("\n");
      }
    }
  }

  return N_interactions;
}

// AmberParm::GetLJparam()
/** Get the LJ 6-12 A and B coefficients (in Amber units) for the given 
  * pair of atoms.
  */
int AmberParm::GetLJparam(double *A, double *B, int atom1, int atom2) {
  int param, index;
  // atype_index = IAC(NATOM)
  // NB_index    = ICO(NTYPES*NTYPES)
  if (LJ_A==NULL || LJ_B==NULL) {
    mprinterr("Error: param file %s does not have LJ A/B coefficients.\n",parmName);
    return 1;
  }
  if (atype_index==NULL || NB_index==NULL) {
    mprinterr("Error: param file %s does not have LJ index information.\n",parmName);
    return 1;
  }
  param = ((ntypes*(atype_index[atom1]-1))+atype_index[atom2])-1; // cpptraj arrays start from 0
  index = NB_index[param]-1;                                      // cpptraj arrays start from 0
  *A = LJ_A[index];
  *B = LJ_B[index];
  return 0;
}

// AmberParm::GetBondParamIdx()
/** Get the bond parameters from the bond_req and bond_rk arrays for
  * the given index.
  */
int AmberParm::GetBondParamIdx(int idxIn, double *rk, double *req) {
  if (bond_rk==NULL || bond_req==NULL) return 1;
  if (idxIn < 0 || idxIn >= numbnd) return 1;
  *rk = bond_rk[idxIn];
  *req = bond_req[idxIn];
  return 0;
}

// AmberParm::GetBondedCutoff()
/** Return bond distance for the two given atoms based on their names.
  */
double AmberParm::GetBondedCutoff(int atom1, int atom2) {
  if (atom1 < 0 || atom1 >= natom) return 0;
  if (atom2 < 0 || atom2 >= natom) return 0;
  return GetBondedCut(names[atom1], names[atom2]);
}

// AmberParm::GetBondParam()
/// Get bond parameters (if they exist) between atom1 and atom2
/** \return 1 if parameters were found, 0 if not.
  */
/*int AmberParm::GetBondParam(double *rk, double *req, int atom1, int atom2) {
  int idx = -2;
  int atom1_3 = atom1 * 3;
  int atom2_3 = atom2 * 3;

  if (bonds==NULL || bondsh==NULL) return 0;
  // Based on name, check bonds or bondsh array, find atom index.
  if (names[atom1][0]=='H' || names[atom2][0]=='H') {
    for (int bond = 0; bond < NbondsWithH * 3; bond += 3) {
      if ( (bondsh[bond] == atom1_3 && bondsh[bond+1] == atom2_3) ||
           (bondsh[bond] == atom2_3 && bondsh[bond+1] == atom1_3)   )
      {
        idx = bondsh[bond+2];
        break;
      }
    }
  } else {
    for (int bond = 0; bond < NbondsWithoutH * 3; bond += 3) {
      if ( (bonds[bond] == atom1_3 && bonds[bond+1] == atom2_3) ||
           (bonds[bond] == atom2_3 && bonds[bond+1] == atom1_3)   )
      {
        idx = bonds[bond+2];
        break;
      }
    }
  }
  // If idx is -2 these atoms are not bonded.
  if (idx==-2) return 0;
  // If idx is -1 atoms are bonded but no parameters for these atoms. Return 
  // default cutoff, no constant.
  if (idx==-1 || bond_rk==NULL || bond_req==NULL) {
    *req = GetBondedCut(names[atom1], names[atom2]);
    *rk = 0.0;
  } else {
    // indices in Amber parms are shifted +1
    idx--;
    *req = bond_req[idx];
    *rk = bond_rk[idx];
  }
  //mprintf("DEBUG:\tAtoms %i and %i are bonded, rk=%lf  req=%lf\n",atom1+1,atom2+1,*rk,*req);
  return 1;
}*/

// AmberParm::SetCharges()
/** Set the atomic charges from the given array. */
int AmberParm::SetCharges(double *chargeIn) {
  if (chargeIn==NULL) return 1;
  if (charge==NULL) charge = new double[ natom ];
  memcpy(charge, chargeIn, natom * sizeof(double));
  return 0;
}

// AmberParm::AmberChargeArray()
/** Charge is in units of electron charge, distance is in angstroms, so 
  * the electrostatic prefactor should be 332. However, since the charges
  * in AmberParm have been converted from Amber charge units, create a new 
  * charged array multiplied by 18.2223. This makes calcs with * Amber-
  * converted charges more accurate at the cost of making calcs with non-Amber 
  * charges less accurate.
  * \param atom_charge Vector to be set with atomic charges * 18.2223
  * \return 0 on success, 1 if no charge information present in parm.
  */
int AmberParm::AmberChargeArray(std::vector<double>& atom_charge) {
  if (charge==NULL) return 1;
  atom_charge.clear();
  atom_charge.resize(natom, ELECTOAMBER);
  for (int atom = 0; atom < natom; atom++)
    atom_charge[atom] *= charge[atom];
  return 0;
}

// AmberParm::AtomCharge()
/** Return charge on given atom. */
double AmberParm::AtomCharge(int atomIn) {
  if (charge==NULL) return 0;
  if (atomIn<0 || atomIn >= natom) return 0;
  return charge[atomIn];
}

// AmberParm::AtomsPerMol()
/** Return number of atoms in given molecule. */
int AmberParm::AtomsPerMol(int molIn) {
  if (atomsPerMol==NULL) return 0;
  if (molIn < 0 || molIn >= molecules) return 0;
  return atomsPerMol[molIn];
}

// ---------- ROUTINES PERTAINING TO SURFACE AREA ------------------------------
// AssignLCPO()
/** Assign parameters for LCPO method. All radii are incremented by 1.4 Ang.
  */
// NOTE: Member function so it can have access to SurfInfo.
void AmberParm::AssignLCPO(SurfInfo *S, double vdwradii, double P1, double P2, 
                           double P3, double P4) {
  S->vdwradii = vdwradii + 1.4;
  S->P1 = P1;
  S->P2 = P2;
  S->P3 = P3;
  S->P4 = P4;
}

// WarnLCPO()
/// Called when the number of bonds to the atom of type atype is not usual.
static void WarnLCPO(char *atype, int atom, int numBonds) {
  mprintf("Warning: Unusual number of bonds for atom %i (%i), type %-2s.\n",
          atom, numBonds, atype);
  mprintf("Using default atom parameters.\n");
}

// AmberParm::SetSurfaceInfo()
/** Set up parameters only used in surface area calcs.
  * LCPO method from:
  *   J. Weiser, P.S. Shenkin, and W.C. Still,
  *   "Approximate atomic surfaces from linear combinations of pairwise
  *   overlaps (LCPO)", J. Comp. Chem. 20:217 (1999).
  * Adapted from gbsa=1 method in SANDER, mdread.f
  * \return the number of solute atoms for which paramters were set. 
  * \return -1 on error.
  */
int AmberParm::SetSurfaceInfo() {
  int *numBonds; // # of bonded neighbors each atom has (LCPO only?)
  int i,atom1,atom2;
  char atype[3];
  atype[2] = '\0';
 
  // If surface info already set up exit 
  if (SurfaceInfo!=NULL) return numSoluteAtoms;
 
  // If no bond information exit
  if (bonds==NULL) {
    mprintf("Error: SetSurfaceInfo(): Parm %s does not contain bond info.\n",parmName);
    return -1;
  } 

  // If no atom type information exit
  if (types==NULL) {
    mprintf("Error: SetSurfaceInfo(): Parm %s does not contain atom type info.\n",
            parmName);
    return -1;
  }
 
  // Get the number of bonded neighbors for each atom
  numBonds = new int[ natom ];
  memset(numBonds, 0, natom * sizeof(int));
  for (i = 0; i < NbondsWithoutH*3; i+=3) {
     atom1 = bonds[i  ] / 3;
     atom2 = bonds[i+1] / 3;
     numBonds[atom1]++;
     numBonds[atom2]++;
  }

  // DEBUG
  //for (i=0; i<natom; i++)
  //  fprintf(stdout,"DEBUG:    Atom %6i_%4s: %2i bonds.\n",i,names[i],numBonds[i]);

  // Only set parameters for solute atoms
  numSoluteAtoms = 0;
  if (firstSolvMol > 0) {
    i = 0;
    while (i < firstSolvMol) numSoluteAtoms += atomsPerMol[i++];
  } else {
    numSoluteAtoms = natom;
  }
  mprintf("[%s] Setting surface parameters for %i solute atoms.\n",parmName,numSoluteAtoms);

  // Set vdw radii and LCPO parameters
  SurfaceInfo = new SurfInfo[ numSoluteAtoms ];
  for (i=0; i < numSoluteAtoms; i++) {
    atype[0] = types[i][0];
    atype[1] = types[i][1];

    if (atype[0]=='C' && atype[1]=='T') {
      switch ( numBonds[i] ) {
        case 1: AssignLCPO(SurfaceInfo+i, 1.70, 0.77887, -0.28063, -0.0012968, 0.00039328); break;
        case 2: AssignLCPO(SurfaceInfo+i, 1.70, 0.56482, -0.19608, -0.0010219, 0.0002658);  break;
        case 3: AssignLCPO(SurfaceInfo+i, 1.70, 0.23348, -0.072627, -0.00020079, 0.00007967); break;
        case 4: AssignLCPO(SurfaceInfo+i, 1.70, 0.00000, 0.00000, 0.00000, 0.00000); break;
        default: WarnLCPO(atype,i,numBonds[i]);
                AssignLCPO(SurfaceInfo+i, 1.70, 0.77887, -0.28063, -0.0012968, 0.00039328);
      }
    } else if (atype[0]=='C' || atype[0]=='c') {
      switch ( numBonds[i] ) {
        case 2: AssignLCPO(SurfaceInfo+i, 1.70, 0.51245, -0.15966, -0.00019781, 0.00016392); break;
        case 3: AssignLCPO(SurfaceInfo+i, 1.70, 0.070344, -0.019015, -0.000022009, 0.000016875); break;
        default: WarnLCPO(atype,i,numBonds[i]);
                AssignLCPO(SurfaceInfo+i, 1.70, 0.77887, -0.28063, -0.0012968, 0.00039328);
      }
    } else if (atype[0]=='O' && atype[1]==' ') {
      AssignLCPO(SurfaceInfo+i, 1.60, 0.68563, -0.1868, -0.00135573, 0.00023743);
    } else if (atype[0]=='O' && atype[1]=='2') {
      AssignLCPO(SurfaceInfo+i, 1.60, 0.88857, -0.33421, -0.0018683, 0.00049372);
    } else if (atype[0]=='O' || atype[0]=='o') {
      switch (numBonds[i]) {
        case 1: AssignLCPO(SurfaceInfo+i, 1.60, 0.77914, -0.25262, -0.0016056, 0.00035071); break;
        case 2: AssignLCPO(SurfaceInfo+i, 1.60, 0.49392, -0.16038, -0.00015512, 0.00016453); break;
        default: WarnLCPO(atype,i,numBonds[i]);
                AssignLCPO(SurfaceInfo+i, 1.60, 0.77914, -0.25262, -0.0016056, 0.00035071);
      }
    } else if (atype[0]=='N' && atype[1]=='3') {
      switch (numBonds[i]) {
        case 1: AssignLCPO(SurfaceInfo+i, 1.65, 0.078602, -0.29198, -0.0006537, 0.00036247); break;
        case 2: AssignLCPO(SurfaceInfo+i, 1.65, 0.22599, -0.036648, -0.0012297, 0.000080038); break;
        case 3: AssignLCPO(SurfaceInfo+i, 1.65, 0.051481, -0.012603, -0.00032006, 0.000024774); break;
        default: WarnLCPO(atype,i,numBonds[i]);
                AssignLCPO(SurfaceInfo+i, 1.65, 0.078602, -0.29198, -0.0006537, 0.00036247);
      }
    } else if (atype[0]=='N' || atype[0]=='n') {
      switch (numBonds[i]) {
        case 1: AssignLCPO(SurfaceInfo+i, 1.65, 0.73511, -0.22116, -0.00089148, 0.0002523); break;
        case 2: AssignLCPO(SurfaceInfo+i, 1.65, 0.41102, -0.12254, -0.000075448, 0.00011804); break;
        case 3: AssignLCPO(SurfaceInfo+i, 1.65, 0.062577, -0.017874, -0.00008312, 0.000019849); break;
        default: WarnLCPO(atype,i,numBonds[i]);
                AssignLCPO(SurfaceInfo+i, 1.65, 0.078602, -0.29198, -0.0006537, 0.00036247);
      }
    } else if (atype[0]=='S' && atype[1]=='H') {
      AssignLCPO(SurfaceInfo+i, 1.90, 0.7722, -0.26393, 0.0010629, 0.0002179);
    } else if (atype[0]=='S' || atype[0]=='s') {
      AssignLCPO(SurfaceInfo+i, 1.90, 0.54581, -0.19477, -0.0012873, 0.00029247);
    } else if (atype[0]=='P' || atype[1]=='p') {
      switch (numBonds[i]) {
        case 3: AssignLCPO(SurfaceInfo+i, 1.90, 0.3865, -0.18249, -0.0036598, 0.0004264); break;
        case 4: AssignLCPO(SurfaceInfo+i, 1.90, 0.03873, -0.0089339, 0.0000083582, 0.0000030381); break;
        default: WarnLCPO(atype,i,numBonds[i]);
          AssignLCPO(SurfaceInfo+i, 1.90, 0.3865, -0.18249, -0.0036598, 0.0004264);
      }
    } else if (atype[0]=='Z') {
      AssignLCPO(SurfaceInfo+i, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000);
    } else if (atype[0]=='H' || atype[0]=='h') {
      AssignLCPO(SurfaceInfo+i, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000);
    } else if (atype[0]=='M' && atype[1]=='G') {
      //  Mg radius = 0.99A: ref. 21 in J. Chem. Phys. 1997, 107, 5422
      //  Mg radius = 1.18A: ref. 30 in J. Chem. Phys. 1997, 107, 5422
      //  Mg radius = 1.45A: Aqvist 1992
      //  The following P1-4 values were taken from O.sp3 with two bonded 
      //  neighbors -> O has the smallest van der Waals radius 
      //  compared to all other elements which had been parametrized
      AssignLCPO(SurfaceInfo+i, 1.18, 0.49392, -0.16038, -0.00015512, 0.00016453);
    } else {
      mprintf("Warning: Using carbon SA parms for unknown atom type %i %2s\n",i,atype);
      AssignLCPO(SurfaceInfo+i, 1.70, 0.51245, -0.15966, -0.00019781, 0.00016392);
    }
  } // END LOOP OVER numSoluteAtoms 

  // DEBUG
  /*
  for (i=0; i<numSoluteAtoms; i++) {
    fprintf(stdout,"%6i %4s: %6.2lf %lf %lf %lf %lf\n",i+1,types[i],SurfaceInfo[i].vdwradii,
    fprintf(stdout,"%6i%6.2lf%12.8lf%12.8lf%12.8lf%12.8lf\n",i+1,SurfaceInfo[i].vdwradii,
            SurfaceInfo[i].P1,SurfaceInfo[i].P2,SurfaceInfo[i].P3,SurfaceInfo[i].P4);
  }
  */
  delete[] numBonds;
  return numSoluteAtoms;
}

// ---------- ROUTINES PERTAINING TO SOLVENT INFO ------------------------------
// AmberParm::IsSolventResname()
/** Return true if the residue name corresponds to solvent. */
bool AmberParm::IsSolventResname(NAME resnameIn) {
  if ( strcmp("WAT ", resnameIn) == 0 ||
       strcmp(" WAT", resnameIn) == 0 ||
       strcmp("HOH ", resnameIn) == 0 ||
       strcmp(" HOH", resnameIn) == 0 ||
       strcmp("TIP3", resnameIn) == 0 
     )
  {
    return true;
  }
  return false;
}

// AmberParm::SetSolventInfo()
/** If atomsPerMol has been read in and firstSolvMol is set, determine solvent 
  * information based on what firstSolvMol is. If firstSolvMol is not set, 
  * determine solvent information by residue name, setting/resetting 
  * atomsPerMol as necessary.
  */
int AmberParm::SetSolventInfo() {
  int molAtom, maskAtom; 

  // Allocate memory
  // Since the number of solvent molecules is not yet known allocate
  // natom for solventMoleculeX arrays. Will be resized after.
  solventMask=new char[ natom ];
  memset(solventMask, 'F', natom * sizeof(char));
  solventMoleculeStart=new int[ natom ];
  solventMoleculeStop=new int[ natom ];
  solventMolecules=0;
  solventAtoms=0;

  // If atomsPerMol is set and firstSolvMol (nspsol) is also set, treat all 
  // the molecules starting with firstSolvMol as solvent.
  if (atomsPerMol!=NULL && firstSolvMol!=-1) {
    molAtom = 0;
    for (int mol=0; mol < molecules; mol++) {
      if (mol+1 >= firstSolvMol) {
        // Add this molecule to the solvent list
        solventAtoms += atomsPerMol[mol];
        for (maskAtom=molAtom; maskAtom < molAtom+atomsPerMol[mol]; maskAtom++)
          solventMask[maskAtom] = 'T';
        solventMoleculeStart[solventMolecules] = molAtom;
        solventMoleculeStop[ solventMolecules] = molAtom+atomsPerMol[mol];
        solventMolecules++;
      }
      molAtom += atomsPerMol[mol];
    }

  // Treat all residues with a recognized solvent name as solvent. This will 
  // reset atomsPerMol from the first solvent molecule on. If atomsPerMol is 
  // not set consider all residues up to the first solvent residue to be in a
  // single molecule.
  } else if (resnums!=NULL) {
    firstSolvMol=-1;
    for (int res=0; res < nres; res++) {
      //mprintf("DEBUG:\tConsidering res %i %4s",res,resnames[res]); 
      if ( IsSolventResname(resnames[res])) {
        // Add this residue to the list of solvent 
        molAtom = resnums[res+1] - resnums[res];
        solventAtoms += molAtom;
        solventMoleculeStart[solventMolecules] = resnums[res];
        solventMoleculeStop[ solventMolecules] = resnums[res+1];
        for (maskAtom=resnums[res]; maskAtom < resnums[res+1]; maskAtom++)
          solventMask[maskAtom] = 'T';
        // If firstSolvMol==-1 this residue is the first solvent molecule 
        if (firstSolvMol==-1) {
          // If atomsPerMol is not yet set up, initialize it. Consider all
          // residues up to this one to be in a single molecule.
          if (atomsPerMol==NULL) {
            // First residue is solvent, all is solvent.
            if (res==0) {
              finalSoluteRes=0;   // Starts from 1, Amber convention
              firstSolvMol=1;     // Starts from 1, Amber convention
              molecules=0;
            } else {
              finalSoluteRes=res; // Starts from 1, Amber convention
              firstSolvMol=2;     // Starts from 1, Amber convention
              molecules=1;
              atomsPerMol = new int[ 1 ];
              atomsPerMol[0] = resnums[res];
            }
          } else { 
            molecules = atomToMolecule(resnums[res]);
            firstSolvMol = molecules + 1; // Starts from 1, Amber convention
          }
        } 
        //mprintf(" solvent mol %i, mol %i\n",solventMolecules,molecules); // DEBUG
        // Update atomsPerMol
        int *tempAtomsPerMol = new int[ molecules + 1];
        if (atomsPerMol!=NULL) {
          memcpy(tempAtomsPerMol, atomsPerMol, molecules * sizeof(int));
          delete[] atomsPerMol;
        }
        atomsPerMol = tempAtomsPerMol;
        atomsPerMol[molecules] = molAtom; 
        solventMolecules++;
        molecules++;
      } else { 
        //mprintf(" not solvent.\n"); // DEBUG
        finalSoluteRes = res + 1; // Starts from 1, Amber convention
      }
    } // End loop over residues
  }

  // DEBUG
  //mprintf("MOLECULE INFORMATION:\n");
  //for (int mol = 0; mol < molecules; mol++)
  //  mprintf("\t%8i %8i\n",mol,atomsPerMol[mol]);

  // Deallocate memory if no solvent 
  if (solventMolecules==0) {
    delete[] solventMask;
    solventMask=NULL;
    delete[] solventMoleculeStart;
    solventMoleculeStart=NULL;
    delete[] solventMoleculeStop;
    solventMoleculeStop=NULL;
    hasSolventInfo = false;

  // Resize the solventMoleculeX arrays
  } else {
    int *tempSMstart = new int[ solventMolecules ];
    memcpy(tempSMstart, solventMoleculeStart, solventMolecules * sizeof(int));
    delete[] solventMoleculeStart;
    solventMoleculeStart = tempSMstart;
    int *tempSMstop = new int[ solventMolecules ];
    memcpy(tempSMstop, solventMoleculeStop, solventMolecules * sizeof(int));
    delete[] solventMoleculeStop;
    solventMoleculeStop = tempSMstop;
    hasSolventInfo = true;
  }

  if (debug>0) {
    mprintf("    %i solvent molecules, %i solvent atoms.\n",
            solventMolecules, solventAtoms);
    if (debug>1)
      mprintf("    FirstSolvMol= %i, FinalSoluteRes= %i\n",firstSolvMol,finalSoluteRes);
  }

  return 0; 
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
// |--------- ROUTINES PERTAINING TO READING PARAMETERS -----------------------|
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
// AmberParm::SetParmName()
void AmberParm::SetParmName(char *parmNameIn, char *parmfileNameIn) {
  // Copy parm filename to parmName. Separate from File.filename in case of stripped parm
  parmName=new char[ strlen(parmNameIn)+1 ];
  strcpy(parmName,parmNameIn);
  parmfileName=new char[ strlen(parmfileNameIn)+1 ];
  strcpy(parmfileName,parmfileNameIn);
}

// AmberParm::CommonSetup()
/** Perform any additional setup required by AmberParm. Should be called once
  * after the AmberParm class is read to in ParmFile. 
  */
// NOTE: Need to ensure this is called ONCE
int AmberParm::CommonSetup(bool bondsearch, bool molsearch) {
  // Create a last dummy residue in resnums that holds natom, which would be
  // the atom number of the next residue if it existed. Atom #s in resnums
  // should correspond with cpptraj atom #s (start from 0) instead of Amber
  // atom #s (start from 1). 
  // Do this to be consistent with PtrajMask selection behavior - saves an 
  // if-then statement.
  int *tempResnums = new int[ nres+1 ];
  memcpy(tempResnums, resnums, nres * sizeof(int));
  delete[] resnums;
  resnums = tempResnums;
  resnums[nres]=natom;
  // DEBUG
  //fprintf(stdout,"==== DEBUG ==== Resnums for %s:\n",parmfile.filename);
  //for (err=0; err<nres; err++) 
  //  fprintf(stdout,"    %i: %i\n",err,resnums[err]);

  // Set up bond information if specified and necessary
  if (bondsearch) {
    if (bonds==NULL && bondsh==NULL && parmCoords!=NULL)
      GetBondsFromCoords(parmCoords);
  }

  // Set up molecule information if specified and necessary
  if (molsearch) {
    if (atomsPerMol==NULL)
      DetermineMolecules();
  }

  // Set up solvent information
  if (SetSolventInfo()) return 1;

  if (debug>0) {
    mprintf("  Number of atoms= %i\n",natom);
    mprintf("  Number of residues= %i\n",nres);
    mprintf("  Number of molecules= %i\n",molecules);
  }

  // Free coords if they were allocated
  if (parmCoords!=NULL) delete[] parmCoords;
  parmCoords=NULL;
  return 0;
}

// ---------- ROUTINES FOR PRINTING PARM INFORMATION ---------------------------
// AmberParm::AtomInfo()
/** Print parm information for atom. */
void AmberParm::AtomInfo(int atom) {
  int res = atomToResidue(atom);
  mprintf("  Atom %i:",atom+1);
  mprintf("[%s]",names[atom]);
  mprintf(" Res %i:",res+1);
  mprintf("[%s]",resnames[res]);
  if (molecules>0)
    mprintf(" Mol %i", atomToMolecule(atom)+1);
  if (types!=NULL)
    mprintf(" Type=[%s]",types[atom]);
  if (charge!=NULL)
    mprintf(" Charge=%lf",charge[atom]);
  if (mass!=NULL)
    mprintf(" Mass=%lf",mass[atom]);
  mprintf("\n");
}

// AmberParm::ParmInfo()
/** Print information about this parm to 1 line. */
void AmberParm::ParmInfo(std::string& ParmTag) {
  if (!ParmTag.empty())
    mprintf(" %i: %s, %i atoms, %i res",pindex,ParmTag.c_str(),natom,nres);
  else if (parmfileName!=NULL)
    mprintf(" %i: %s, %i atoms, %i res",pindex,parmfileName,natom,nres);
  else
    mprintf(" %i: %s, %i atoms, %i res",pindex,parmName,natom,nres);
  if (boxType==NOBOX)
    mprintf(", no box");
  else if (boxType==ORTHO)
    mprintf(", ortho. box");
  else if (boxType==NONORTHO)
    mprintf(", non-ortho. box");
  if (molecules>0)
    mprintf(", %i mol",molecules);
  if (solventMolecules>0)
    mprintf(", %i solvent mol",solventMolecules);
  if (parmFrames>0)
    mprintf(", %i frames",parmFrames);
  mprintf("\n");
}

// AmberParm::Summary()
/** Print a summary of atoms, residues, molecules, and solvent molecules
  * in this parm.
  */
void AmberParm::Summary() {
  mprintf("              Topology contains %i atoms.\n",this->natom);
  mprintf("                                %i residues.\n",this->nres);
  int number_of_bonds = NbondsWithH + NbondsWithoutH;
  mprintf("                                %i bonds.\n",number_of_bonds);
  if (this->molecules>0)
    mprintf("                                %i molecules.\n",this->molecules);
  if (this->solventMolecules>0) {
    mprintf("                                %i solvent molecules.\n",this->solventMolecules);
    mprintf("                  First solvent molecule is %i\n",this->firstSolvMol);
    mprintf("                  Final solute residue is %i\n",this->finalSoluteRes);
  }
}

// AmberParm::PrintBondInfo()
/** Print information contained in bonds and bondsh arrays. */
void AmberParm::PrintBondInfo() {
  int atom1,atom2,atomi;
  if ((NbondsWithH + NbondsWithoutH) <= 0) {
    mprintf("NO BOND INFORMATION IN PRMTOP\n");
    return;
  }
  if (NbondsWithH>0) {
    mprintf("%i BONDS TO HYDROGEN:\n",NbondsWithH);
    for (int ibond=0; ibond < NbondsWithH * 3; ibond += 3) {
      atom1 = (bondsh[ibond  ]/3) + 1;
      atom2 = (bondsh[ibond+1]/3) + 1;
      atomi = (bondsh[ibond+2]  );
      mprintf("\tAtom %i to %i, %i\n",atom1,atom2,atomi);
    }
  }
  if (NbondsWithoutH>0) {
    mprintf("%i BONDS TO NON-HYDROGEN:\n",NbondsWithoutH);
    for (int ibond=0; ibond < NbondsWithoutH * 3; ibond += 3) {
      atom1 = (bonds[ibond  ]/3) + 1;
      atom2 = (bonds[ibond+1]/3) + 1;
      atomi = (bonds[ibond+2]  );
      mprintf("\tAtom %i to %i, %i\n",atom1,atom2,atomi);
    }
  }
}

// AmberParm::PrintMoleculeInfo()
/** Print information on molecules. */
void AmberParm::PrintMoleculeInfo() {
  int atomcount = 0;
  int resid;
  char rtemp[32];
  if (molecules==0 || atomsPerMol==NULL) {
    mprintf("NO MOLECULE INFORMATION IN PRMTOP\n");
    return;
  }
  mprintf("MOLECULES:\n");
  for (int mol=0; mol < molecules; mol++) {
    resid = atomToResidue(atomcount);
    ResName(rtemp,resid);
    mprintf("\tMolecule %i, %i atoms, first residue %s\n",mol+1,atomsPerMol[mol],rtemp);
    atomcount += atomsPerMol[mol];
  }
}

// AmberParm::PrintResidueInfo()
/** Print information on residues. */
void AmberParm::PrintResidueInfo() {
  mprintf("RESIDUES:\n");
  for (int res = 0; res < nres; res++) 
    mprintf("\tResidue %i, %s, first atom %i\n",res+1,resnames[res],resnums[res]+1);
}

// ---------- ROUTINES FOR CONVERTING ATOM/RESIDUE/MOLCULE NUMBERS -------------
// AmberParm::atomToResidue()
/** Given an atom number, return corresponding residue number. */
int AmberParm::atomToResidue(int atom) {
  int i;
  if (atom >= 0 && atom < natom) {
    for (i = 0; i < nres; i++)
      if ( atom>=resnums[i] && atom<resnums[i+1] )
        return i;
  }
  return -1;
}

// AmberParm::atomToMolecule()
/** Given an atom number, return corresponding molecule number. */
int AmberParm::atomToMolecule(int atom) {
  int a = 0;
  if (atom >= 0 && atom < natom) {
    for (int i = 0; i < molecules; i++) {
      a += atomsPerMol[i];
      if (atom < a) return i;
    }
  }
  return -1;
}

// AmberParm::atomToSolventMolecule()
/** Given an atom number, return corresponding solvent molecule. */
// NOTE: Could this be achieved with atomToMolecule and solventMask?
int AmberParm::atomToSolventMolecule(int atom) {
  int i, atom1;
  if (atom >= 0 && atom < natom) {
    atom1 = atom + 1; 
    for (i = 0; i < molecules; i++) {
      if (atom1 <= solventMoleculeStart[i])
        return -1;
      else if (atom1>solventMoleculeStart[i] && atom1<=solventMoleculeStop[i])
        return i;
    }
  }
  return -1;
}

// ---------- ROUTINES PERTAINING TO BOND INFORMATION --------------------------
// AmberParm::ResetBondInfo()
/** Reset the bonds and bondsh arrays, as well as NBONH and NBONA */
void AmberParm::ResetBondInfo() {
  if (bonds!=NULL) delete[] bonds;
  bonds=NULL;
  if (bondsh!=NULL) delete[] bondsh;
  bondsh=NULL;
  NbondsWithH=0;
  NbondsWithoutH=0;
  bondInfo.Reset();
}

// AmberParm::AddBond()
/** Add bond info for the two atoms. Attempt to identify if it is a bond
  * to hydrogen or not based on names. The atom numbers should start from 0.
  * Atom indices in bond arrays are * 3.
  */ 
int AmberParm::AddBond(int atom1, int atom2, int icb) {
  bool isH=false;
  int bondidx;
  if (atom1<0 || atom2<0 || atom1>=natom || atom2>=natom) return 1;
  if (names!=NULL) {
    if (names[atom1][0]=='H') isH=true;
    if (names[atom2][0]=='H') isH=true;
  }
  if (isH) {
    //NbondsWithH = values[NBONH];
    bondidx = NbondsWithH * 3;
    int *tempbondsh = new int[ bondidx+3 ];
    memcpy(tempbondsh, bondsh, bondidx * sizeof(int));
    delete[] bondsh;
    bondsh = tempbondsh;
    bondsh[bondidx  ] = atom1 * 3;
    bondsh[bondidx+1] = atom2 * 3;
    bondsh[bondidx+2] = icb;
    NbondsWithH++;
  } else {
    //NbondsWithoutH = values[NBONA];
    bondidx = NbondsWithoutH * 3;
    int *tempbonds = new int[ bondidx+3 ];
    memcpy(tempbonds, bonds, bondidx * sizeof(int));
    delete[] bonds;
    bonds = tempbonds;
    bonds[bondidx  ] = atom1 * 3;
    bonds[bondidx+1] = atom2 * 3;
    bonds[bondidx+2] = icb;
    NbondsWithoutH++;
  }
  return 0;
}

// AmberParm::GetBondsFromCoords()
/** Given an array of coordinates X0Y0Z0X1Y1Z1...XNYNZN determine which
  * atoms are bonded via distance search. First check for bonds within
  * residues, then check for bonds between adjacent residues. Adjacent
  * residues in different molecules are not considered.
  */
// NOTE: Speedup by precalc elements, use grid, get cutoff^2
void AmberParm::GetBondsFromCoords(double *inputCoords) {
  int res, startatom, stopatom, midatom,atom1, atom2, idx1, idx2, stopResnum;
  double D2, cutoff2;
  int *resmols;
  int firstSolventResidue = -1;

  // NOTE: Can either attempt to be very accurate by determining atomic
  // element from atom name and getting the cutoff that way, use input radii
  // to determine an appropriate cutoff (the largest radius scaled by some
  // factor), or just use a fixed cutoff (which is very fast). Use cutoff^2
  // to compare to dist^2 and avoid the sqrt op.
  //cutoff2 = 0.99920016; // (0.833 * 1.2)^2

  if (inputCoords==NULL) return;
  mprintf("\t%s: determining bond info from distances.\n",parmName);
  // Determine bonds within residues.
  for (res = 0; res < nres; res++) {
    // Check if this residue is the first solvent molecule.
    if (firstSolventResidue==-1 && IsSolventResname(resnames[res]))
      firstSolventResidue = res;
    startatom = resnums[res];
    stopatom = resnums[res+1];
    //mprintf("\t\tDetermining bonds within residue %i\n",res);
    for (atom1 = startatom; atom1 < stopatom - 1; atom1++) {
      idx1 = atom1 * 3;
      for (atom2 = atom1 + 1; atom2 < stopatom; atom2++) {
        idx2 = atom2 * 3;
        D2 = DIST2_NoImage(inputCoords + idx1, inputCoords + idx2);
        //mprintf("\t\tGetting cutoff for [%s] - [%s]\n",names[atom1],names[atom2]);
        cutoff2 = GetBondedCut(names[atom1],names[atom2]);
        cutoff2 *= cutoff2; // Op '*' less expensive than sqrt
//        if (debug>0) {
//          if (debug==1) {
//            if (D2<cutoff2) mprintf("\tBOND: %s %i to %s %i\n",names[atom1],atom1+1,
//                                    names[atom2],atom2+1);
//          } else if (debug > 1) {
//            mprintf("Distance between %s %i and %s %i is %lf, cutoff2 %lf",names[atom1],atom1+1,
//                    names[atom2],atom2+1,D2,cutoff2);
//            if (D2<cutoff2) mprintf(" bonded!");
//            mprintf("\n");
//          }
//        }
        if (D2 < cutoff2) {
          AddBond(atom1,atom2,-1);
          // Test: Once a bond has been made to hydrogen move on
          if (names[atom1][0]=='H') break;
        }
      }
    }
  }

  // Dont want to check for bonds between residues once we are in a 
  // solvent region, so set the last residue to be searched to the final
  // solute residue.
  if (firstSolventResidue!=-1)
    stopResnum = firstSolventResidue;
  else
    stopResnum = nres;

  // If atomsPerMol has been set up, create an array that will contain the 
  // molecule number of each residue.
  resmols = new int[ stopResnum ];
  if (atomsPerMol!=NULL) {
    int molnum = 0;
    int atotal = atomsPerMol[0];
    for (res = 0; res < stopResnum; res++) {
      resmols[res] = molnum;
      if (resnums[res+1] >= atotal) {
        molnum++;
        if (molnum >= molecules) break;
        atotal += atomsPerMol[molnum];
      }
    }
  } else
    memset(resmols, 0, stopResnum * sizeof(int));
  // DEBUG
  //for (res = 0; res < stopResnum; res++) 
  //  mprintf("DEBUG\tRes %8i %4s Mol %8i\n",res+1,resnames[res],resmols[res]);

  // Determine bonds between adjacent residues. 
  for (res = 1; res < stopResnum; res++) {
    // Dont check for bonds between residues that are in different molecules
    if (resmols[res-1] != resmols[res]) continue;
    startatom = resnums[res-1];
    midatom = resnums[res];
    stopatom = resnums[res+1];
    //mprintf("\t\tDetermining bonds between residues %i and %i\n",res-1,res);
    for (atom1 = startatom; atom1 < midatom; atom1++) {
      idx1 = atom1 * 3;
      // Test: Hydrogen should not bond across residues
      if (names[atom1][0]=='H') continue;
      for (atom2 = midatom; atom2 < stopatom; atom2++) {
        idx2 = atom2 * 3;
        // Test: Hydrogen should not bond across residues
        if (names[atom2][0]=='H') continue;
        D2 = DIST2_NoImage(inputCoords + idx1, inputCoords + idx2);
        cutoff2 = GetBondedCut(names[atom1],names[atom2]);
        cutoff2 *= cutoff2;
        if (D2 < cutoff2) AddBond(atom1,atom2,-1);
      } 
    }
  }
  
  delete[] resmols;

  mprintf("\t%s: %i bonds to hydrogen, %i other bonds.\n",parmName,NbondsWithH,NbondsWithoutH);
}

// AmberParm::BondArray()
/// Set the given array with bond information.
/** Set the given vector up with information from the bonds and bondsh arrays.
  * The vector will have format Bond1Atom1, Bond1Atom2, Bond2Atom1, Bond2Atom2
  * etc. First set bond without H, then bonds with H. Currently only used by
  * Traj_Mol2file?
  * \param bondarray vector to be set up with bond information
  * \return 0 if successfully set up, 1 if error occurs.
  */
int AmberParm::BondArray(std::vector<int> &bondarray) {
  if (bonds==NULL && bondsh==NULL) return 1;
  bondarray.clear();
  bondarray.reserve( (NbondsWithoutH*2) + (NbondsWithH*2) );
  // NOTE: Atom #s in the bonds and bondh arrays are * 3
  if (bonds!=NULL) {
    int nb = NbondsWithoutH * 3;
    for (int bond = 0; bond < nb; bond+=3) {
      bondarray.push_back(bonds[bond]/3);
      bondarray.push_back(bonds[bond+1]/3);
    }
  }
  if (bondsh!=NULL) {
    int nbh = NbondsWithH * 3;
    for (int bond = 0; bond < nbh; bond+=3) {
      bondarray.push_back(bondsh[bond]/3);
      bondarray.push_back(bondsh[bond+1]/3);
    }
  }
  return 0;
}

// AmberParm::BondArrayWithParmIdx()
/// Set the given array with bond information including Amber parm bond index.
/** Set the given vector up with information from the bonds and bondsh arrays.
  * The vector will have format Bond1Atom1, Bond1Atom2, Bond1Idx, Bond2Atom1, 
  * Bond2Atom2, Bond2Idx, etc. Ensure that BondXAtom1 is always < BondXAtom2
  * for easy use in pairwise loops. First set bond without H, then bonds with H.
  * Currently only used by Action_CheckStructure?
  * \param bondarray vector to be set up with bond information
  * \return 0 if successfully set up, 1 if error occurs.
  */
int AmberParm::BondArrayWithParmIdx(std::vector<int> &bondarray) {
  int batom1, batom2, tempatom;
  if (bonds==NULL && bondsh==NULL) return 1;
  bondarray.clear();
  int nb = NbondsWithoutH * 3;
  int nbh = NbondsWithH * 3;
  bondarray.reserve( nb + nbh );
  // NOTE: Atom #s in the bonds and bondh arrays are * 3
  // NOTE: Indices initially start from 1, so subtract 1
  if (bonds!=NULL) {
    for (int bond = 0; bond < nb; bond+=3) {
      batom1 = bonds[bond]/3;
      batom2 = bonds[bond+1]/3;
      if (batom2 < batom1) {
        tempatom = batom1;
        batom1 = batom2;
        batom2 = tempatom;
      }
      bondarray.push_back(batom1);
      bondarray.push_back(batom2);
      bondarray.push_back(bonds[bond+2]-1);
    }
  }
  if (bondsh!=NULL) {
    for (int bond = 0; bond < nbh; bond+=3) {
      batom1 = bondsh[bond]/3;
      batom2 = bondsh[bond+1]/3;
      if (batom2 < batom1) {
        tempatom = batom1;
        batom1 = batom2;
        batom2 = tempatom;
      }
      bondarray.push_back(batom1);
      bondarray.push_back(batom2);
      bondarray.push_back(bondsh[bond+2]-1);
    }
  }
  return 0;
}

// ---------- ROUTINES PERTAINING TO BONDS IN BONDINFO CLASS -------------------
// AmberParm::SetupBondInfo()
/// Set up the BondInfo structure from the bonds and bondsh arrays.
// NOTE: This can eventually be used to replace the bonds and bondsh arrays,
// which are not optimal for things like determining whether one atom is
// bonded to another, or getting a list of atoms bonded to a certain atom.
int AmberParm::SetupBondInfo() {
  // If already set up exit now.
  if (bondInfo.HasBeenSetup()) return 0;
  if (bonds==NULL && bondsh==NULL) {
    mprinterr("Error: SetupBondInfo: No bond information present.\n");
    return 1;
  }
  mprintf("\t%s: Setting up %i bonds.\n",parmName,NbondsWithH+NbondsWithoutH);

  bondInfo.Setup(natom);

  // Set max valences
  bondInfo.SetValences(names);

  // Go through the bonds and bondsh arrays
  bondInfo.SetBondsFromAmberArray(bondsh, NbondsWithH);
  bondInfo.SetBondsFromAmberArray(bonds,  NbondsWithoutH);
  return 0;
}

// AmberParm::GetBondedAtomIdx()
/// If given atom name matches one bonded to given atom#, return its atom# 
/** This routine assumes the BondInfo structure has been set up from a call
  * to SetupBondInfo. Use this to determine whether the given atom name
  * matches one of the atoms bonded to the given atom number.
  * \return the atom# of bonded atom, -1 if not found.
  */
int AmberParm::GetBondedAtomIdx(int atomIn, const char *bondedAtomName) {
  //int bondList[MAXNUMBONDS]; // Defined in Bonds.h
  std::vector<int> bondList;

  bondInfo.GetListOfBondedAtoms(atomIn, bondList);
  for (std::vector<int>::iterator bondedAtom = bondList.begin();
                                  bondedAtom != bondList.end();
                                  bondedAtom++) {
    if (strcmp(names[ *bondedAtom ], bondedAtomName)==0) return *bondedAtom;
  }
  return -1;
}

// AmberParm::GetBondedHatoms()
/// Set up a vector containing the atom #s of any H atoms bonded to this atom
/** \return 1 if there are bonded H atoms, 0 if no bonded H atoms.
  */
int AmberParm::GetBondedHatoms(int atomIn, std::vector<int>& HatomList) {
  std::vector<int> bondList;

  HatomList.clear();  
  bondInfo.GetListOfBondedAtoms(atomIn, bondList);
  for (std::vector<int>::iterator bondedAtom = bondList.begin();
                                  bondedAtom != bondList.end();
                                  bondedAtom++) {
    // DEBUG
    //mprintf("\t\tAtom %i@%s bonded to %i@%4s",atomIn+1,names[atomIn],
    //        *bondedAtom+1,names[*bondedAtom]);
    if (names[ *bondedAtom ][0] == 'H') {
      HatomList.push_back( *bondedAtom );
      //mprintf(" HYDROGEN!");
    }
    //mprintf("\n");
  }
  if (HatomList.empty()) return 0;
  return 1; 
}

// AmberParm::MaskOfAtomsAroundBond()
/// Just a wrapper for the same routine in BondInfo
/** Requires the internal BondInfo structure to be set up. */
int AmberParm::MaskOfAtomsAroundBond(int atom1, int atom2, std::vector<char> &Selected) {
  return bondInfo.MaskOfAtomsAroundBond(atom1,atom2,Selected);
}

// AmberParm::DetermineMolecules()
/** Given that bonding information for the parm has been set up, attempt
  * to determine how many molecules (i.e. entities that are not covalently
  * bonded) there are.
  */
int AmberParm::DetermineMolecules() {
  if (SetupBondInfo()) {
    mprinterr("Error: DetermineMolecules: No bond information set up.\n");
    return 1;
  }
  mprintf("\t%s: Determining molecule information from bonds.\n",parmName);

  // Determine molecules
  atomsPerMol = bondInfo.DetermineMolecules(&molecules);
 
  //mol.PrintBonds();

  return 0;
}

// AmberParm::SetupExcludedAtoms()
/** Setup excluded atoms list and numex based on bond info */
int AmberParm::SetupExcludedAtoms() {
  if (SetupBondInfo()) {
    mprinterr("Error: SetupExcludedAtoms: No bond information in parm %s.\n",parmName);
    return 1;
  }
  if (numex!=NULL) delete[] numex;
  numex = new int[ natom ];
  if (excludedAtoms!=NULL) delete[] excludedAtoms;
  excludedAtoms = bondInfo.DetermineExcludedAtoms(numex, &nnb);
  return 0;
}

// ---------- ROUTINES PERTAINING TO MODIFYING THE PARM ------------------------
// SetupSequentialArray()
/// Given an atom map and old sequential array, set up new sequential array.
/** Given an array with format [I0J0...X0][I1J1...] where the entries up to
  * X are atom# * 3 and entry X is an index, and an atom map array
  * with format Map[oldAtom]=newAtom, create a new array that contains
  * only entries for which all atoms are present. Can be used for the
  * bond, angle, and dihedral arrays.
  */
// NOTE: Set up atom map to be atom*3??
static int *SetupSequentialArray(int *atomMap, int oldN, int Nsequence, 
                                 int *oldArray, int *newN) {
  int *newArray = NULL;
  std::vector<int> newatoms;
  std::vector<int> newsign;
  
  if (atomMap==NULL || oldArray==NULL) return NULL;
  // Actual # entries in oldArray
  int oldNX = oldN * Nsequence;
  // # of coord indices in each entry; last entry is a parameter index
  int Nsequence1 = Nsequence - 1;
  // Set initial size of new array to that of old array
  newArray = new int[ oldN * Nsequence ];
  int newNX = 0;
  // NOTE: Make newatoms static?
  newatoms.resize(Nsequence, 0);
  newsign.resize(Nsequence1, 1);
  // Go through old array, use atomMap to determine what goes into newParm
  for (int oldi=0; oldi < oldNX; oldi += Nsequence) {
    // Check that atoms 0 to Nsequence exist in newParm. If any of the atoms
    // do not exist in newParm bail.
    int newatm = -1;
    bool reverseOrder = false;
    for (int sequencei = 0; sequencei < Nsequence1; sequencei++) {
      int arrayAtom = oldArray[oldi+sequencei];
      // For dihedrals the atom # can be negative. Convert to positive
      // for use in the atom map.
      if (arrayAtom < 0) {
        newatm = atomMap[ -arrayAtom / 3 ];
        newsign[sequencei] = -1;
        // For improper/multi-term dihedrals atom index 0 cannot be the 3rd
        // or 4th position since there is no such thing as -0.
        if (newatm == 0 && (sequencei == 2 || sequencei == 3))
          reverseOrder = true;
      } else {
        newatm = atomMap[  arrayAtom / 3 ];
        newsign[sequencei] = 1;
      }
      // New atom # of -1 means atom was removed - exit loop now.
      if (newatm == -1) break;
      // Atom exists - store.
      newatoms[sequencei] = newatm;
    }
    // If newatm is -1 here that means it didnt exist in newParm for this
    // sequence. Skip the entire sequence.
    if (newatm==-1) continue;
    // Store the final number of the sequence, which is an index
    newatoms[Nsequence1] = oldArray[oldi+Nsequence1];
    // Place the atoms in newatoms in newArray
    if (!reverseOrder) {
      for (int sequencei = 0; sequencei < Nsequence1; sequencei++) 
        newArray[newNX+sequencei] = newatoms[sequencei] * 3 * newsign[sequencei];
    } else {
      int ridx = Nsequence1 - 1;
      for (int sequencei = 0; sequencei < Nsequence1; sequencei++)
        newArray[newNX+sequencei] = newatoms[ridx--] * 3 * newsign[sequencei];
    }
    // Place the index in newatoms in newArray
    newArray[newNX+Nsequence1] = newatoms[Nsequence1];
    // Increment new counter
    newNX += Nsequence;
  }
  // Resize newArray
  int *temparray = new int[ newNX ];
  memcpy(temparray, newArray, newNX * sizeof(int));
  delete[] newArray;
  newArray = temparray;
  
  *newN = newNX / Nsequence;
  return newArray;
}

// SetupBondArray()
/// Given an atom map and new parm, set up bond array
// NOTE: Set up atom map to be atom*3??
// NOTE: May be obsolete 
static int *SetupBondArray(int *atomMap, int oldN3, int *oldBonds, int *newN) {
  int *bonds;
  int N3, i, atom1, atom2;

  if (atomMap==NULL || oldBonds==NULL) return NULL;
  bonds=NULL;
  N3=0;
  // Go through Bonds with/without H, use atomMap to determine what goes into newParm
  for (i=0; i < oldN3; i+=3) {
    // Check that atom1 and atom2 exist in newParm
    // In the bond arrays atom nums are multiplied by 3
    atom1 = atomMap[ oldBonds[i]/3   ];
    atom2 = atomMap[ oldBonds[i+1]/3 ];
    if ( atom1!=-1 && atom2!=-1 ) {
      // Put new atom 1 and new atom 2 in newParm array
      int *tempbonds = new int[ N3+3 ];
      memcpy(tempbonds, bonds, N3 * sizeof(int));
      delete[] bonds;
      bonds = tempbonds; 
      bonds[N3]   = atom1 * 3;
      bonds[N3+1] = atom2 * 3;
      bonds[N3+2] = oldBonds[i+2];
      N3+=3;
    }
  }
  
  *newN = N3 / 3;
  return bonds;
}

// AmberParm::modifyStateByMap()
/** Currently only intended for use with AtomMap.
  * This routine will create a new amber parm (newParm) base on the
  * current amber parm (this), mapping atoms in newParm to atoms
  * in this based on the given atom map.
  * NOTE: There is no guarantee that atoms that were contiguous in 
  *       this parm will be contiguous in the old parm since this is not
  *       currently enforced by AtomMap; therefore the residue information
  *       will probably be shot unless there is only 1 residue. 
  * NOTE: Molecule, solvent info etc is not copied over.
  */
AmberParm *AmberParm::modifyStateByMap(int *AMap) {
  AmberParm *newParm;
  int j=0;
  int *ReverseMap;

  newParm = new AmberParm();
  newParm->SetDebug(debug);
  // Allocate space for arrays and perform initialization
  newParm->names    = new NAME[ natom ];
  if (this->types!=NULL)
    newParm->types    = new NAME[ natom ];
  if (this->charge!=NULL)
    newParm->charge   = new double[ natom ];
  if (this->mass!=NULL)
    newParm->mass     = new double[ natom ];
  newParm->resnames = new NAME[ nres ];
  newParm->resnums  = new int[ nres+1 ];
  // Need reverse of AMap, Map[tgt atom] = ref atom for setting up bonds
  ReverseMap = new int[ natom ];

  // Loop over all atoms in this parm, map them to new parm
  for (int i=0; i < this->natom; i++) {
    j = AMap[i];
    ReverseMap[j] = i;
    strcpy(newParm->names[i], this->names[j]);
    if (this->types!=NULL)  strcpy(newParm->types[i], this->types[j]);
    if (this->charge!=NULL) newParm->charge[i] =      this->charge[j];
    if (this->mass!=NULL)   newParm->mass[i]   =      this->mass[j];
  }

  // Copy residue info. If > 1 residue the copy will likely not be correct.
  if (this->nres>1) {
    mprintf("WARNING: modifyStateByMap: %s has > 1 residue, modified parm residue info\n",parmName);
    mprintf("         will most likely not be correct!\n");
  }
  for (int res=0; res<this->nres; res++) {
    strcpy(newParm->resnames[res],this->resnames[res]);
    newParm->resnums[res] = this->resnums[res];
  }
  // Fix up IPRES
  newParm->resnums[this->nres] = this->natom;

  // Set up bond arrays
  newParm->bondsh = SetupBondArray(ReverseMap, this->NbondsWithH*3, this->bondsh,
                                   &(newParm->NbondsWithH));
  newParm->bonds  = SetupBondArray(ReverseMap, this->NbondsWithoutH*3, this->bonds,
                                   &(newParm->NbondsWithoutH));
  // Clear reverse map
  delete[] ReverseMap; 

  // Set up new parm information
  newParm->natom = this->natom;
  newParm->nres = this->nres;
  newParm->parmFrames = this->parmFrames;

  // Give mapped parm the same pindex as original parm
  newParm->pindex = this->pindex;

  // Copy box information
  for (int i=0; i<6; i++)
    newParm->Box[i] = this->Box[i];
  newParm->boxType=this->boxType;

  return newParm;
}

// AmberParm::modifyStateByMask()
/**  The goal of this routine is to create a new AmberParm (newParm)
  *  based on the current AmberParm (this), deleting atoms that are
  *  not in the Selected array.
  */
// NOTE: Make all solvent/box related info dependent on IFBOX only?
// NOTE: Eventually convert so atom mask is passed in.
AmberParm *AmberParm::modifyStateByMask(std::vector<int> &Selected, char *prefix) {
  AmberParm *newParm;
  int selected;
  int i, ires, imol; 
  int j, jres, jmol;
  int curres, curmol; 
  int *atomMap; // Convert atom # in oldParm to newParm; -1 if atom is not in newParm
  int Nselected = (int) Selected.size();

  // Allocate space for the new state
  newParm = new AmberParm(); 
  newParm->SetDebug(debug);

  // Allocate space for arrays and perform initialization
  atomMap = new int[ natom ];
  for (i=0; i<this->natom; i++) atomMap[i]=-1;
  newParm->names = new NAME[ Nselected ];
  if (this->types!=NULL)
    newParm->types = new NAME[ Nselected ];
  if (this->charge!=NULL)
    newParm->charge = new double[ Nselected ];
  if (this->mass!=NULL)
    newParm->mass = new double[ Nselected ];
  if (this->atype_index!=NULL)
    newParm->atype_index = new int[ Nselected ];
  if (this->at_num!=NULL)
    newParm->at_num = new int[ Nselected ];
  if (this->gb_radii!=NULL)
    newParm->gb_radii = new double[ Nselected ];
  if (this->gb_screen!=NULL) 
    newParm->gb_screen = new double[ Nselected ];
  // The following arrays could eventually be determined somehow, however
  // at the moment its not worth it since sander doesnt use the info
  // anymore anyway.
  if (this->itree!=NULL)
    newParm->itree = new NAME[ Nselected ];
  if (this->join_array!=NULL)
    newParm->join_array = new int[ Nselected ];
  if (this->irotat!=NULL)
    newParm->irotat = new int[ Nselected ];
  // Arrays for which final size will not be known
  newParm->resnames = new NAME[ nres ];
  newParm->resnums  = new int[ nres+1 ];
  // Other arrays
  if (this->radius_set!=NULL) {
    newParm->radius_set = new char[ strlen(this->radius_set)+1 ];
    strcpy(newParm->radius_set, this->radius_set);
  } 
  if (this->molecules>0) 
    newParm->atomsPerMol = new int[ molecules ];

  // Set stripped parm name based on prefix: <prefix>.<oldparmname>
  // If no prefix given set name as: strip.<oldparmname>
  if (prefix!=NULL) {
    newParm->parmName=new char[ strlen(this->parmName)+strlen(prefix)+2 ];
    sprintf(newParm->parmName,"%s.%s",prefix,this->parmName);
  } else {
    newParm->parmName=new char[ strlen(this->parmName)+7];
    sprintf(newParm->parmName,"strip.%s",this->parmName);
  }

  // Set first solvent molecule to -1 for now. If there are no solvent 
  // molecules left in newParm after strip it will be set to 0.
  newParm->firstSolvMol=-1;

  j = 0; 
  jres = -1; jmol = -1;
  ires = -1; imol = -1;

  // Loop over Selected atoms and set up information for the newstate if the atom is 
  // not to be deleted...
  for (selected=0; selected < Nselected; selected++) {
    // i = old atom #, j = new atom number
    i = Selected[selected];          // Atom to be kept from oldParm
    curres = this->atomToResidue(i); // Residue number of atom in oldParm
    atomMap[i]=j;                    // Store this atom in the atom map
    // Copy over atom information
    strcpy(newParm->names[j], this->names[i]);
    if (this->types!=NULL)       strcpy(newParm->types[j], this->types[i]);
    if (this->charge!=NULL)      newParm->charge[j]      = this->charge[i];
    if (this->mass!=NULL)        newParm->mass[j]        = this->mass[i];
    if (this->atype_index!=NULL) newParm->atype_index[j] = this->atype_index[i];
    if (this->at_num!=NULL)      newParm->at_num[j]      = this->at_num[i];
    if (this->gb_radii!=NULL)    newParm->gb_radii[j]    = this->gb_radii[i];
    if (this->gb_screen!=NULL)   newParm->gb_screen[j]   = this->gb_screen[i];
    if (this->itree!=NULL)       strcpy(newParm->itree[j], this->itree[i]);
    if (this->join_array!=NULL)  newParm->join_array[j]  = this->join_array[i];
    if (this->irotat!=NULL)      newParm->irotat[j]      = this->irotat[i];

    // Check to see if we are in the same residue or not and copy relevant information
    if (ires == -1 || ires != curres) {
      jres++;
      strcpy(newParm->resnames[jres], this->resnames[curres]);
      newParm->resnums[jres] = j;
      ires = curres;
    }

    // Check to see if we are in the same molecule or not and increment #atoms in molecule
    if (this->molecules>0) {
      curmol = this->atomToMolecule(i);
      if (imol == -1 || imol != curmol) {
        jmol++;
        newParm->atomsPerMol[jmol]=1;
        imol = curmol;
      } else {
        newParm->atomsPerMol[jmol]++;
      }
    }

    // If we are keeping this atom and it belongs to a solvent molecule and 
    // the first solvent atom has not been set, set it.
    if (this->solventMolecules>0 && this->solventMask[i]=='T' && newParm->firstSolvMol<0) {
      newParm->firstSolvMol = jmol + 1;
      newParm->finalSoluteRes = jres;
    }

    // Increment the new atom counter
    j++;
  } // End loop over oldParm Selected atoms 

  // Set up bond arrays
  newParm->bondsh = SetupSequentialArray(atomMap, this->NbondsWithH, 3, 
                                         this->bondsh, &(newParm->NbondsWithH));
  newParm->bonds  = SetupSequentialArray(atomMap, this->NbondsWithoutH, 3, 
                                         this->bonds,  &(newParm->NbondsWithoutH));
  // Set up angle arrays
  newParm->anglesh = SetupSequentialArray(atomMap, this->NanglesWithH, 4,
                                          this->anglesh, &(newParm->NanglesWithH));
  newParm->angles  = SetupSequentialArray(atomMap, this->NanglesWithoutH, 4,
                                          this->angles,  &(newParm->NanglesWithoutH));
  // Set up dihedral arrays
  newParm->dihedralsh = SetupSequentialArray(atomMap, this->NdihedralsWithH, 5,
                                             this->dihedralsh, &(newParm->NdihedralsWithH));
  newParm->dihedrals  = SetupSequentialArray(atomMap, this->NdihedralsWithoutH, 5,
                                             this->dihedrals,  &(newParm->NdihedralsWithoutH));
  delete[] atomMap;

  // NOTE: Since in the above arrays the indices have survived intact we can 
  //       just include direct copies of all the constants arrays for now. 
  //       May want to cull this later.
  // Nonbonded parm index; since this depends on the atom type index and
  // those entries were not changed this is still valid
  if (NB_index!=NULL) {
    newParm->ntypes = this->ntypes;
    newParm->NB_index = new int[ ntypes * ntypes ];
    memcpy(newParm->NB_index, this->NB_index, ntypes*ntypes * sizeof(int));
    // Lennard-jones params
    int ljmax = ntypes*(ntypes+1)/2;
    if (this->LJ_A!=NULL && this->LJ_B!=NULL) {
      newParm->LJ_A = new double[ ljmax ];
      memcpy(newParm->LJ_A, this->LJ_A, ljmax * sizeof(double));
      newParm->LJ_B = new double[ ljmax ];
      memcpy(newParm->LJ_B, this->LJ_B, ljmax * sizeof(double));
    }
    // LJ 10-12 params
    newParm->nphb = this->nphb;
    if (this->nphb>0) {
      newParm->asol = new double[ nphb ];
      memcpy(newParm->asol, this->asol, nphb * sizeof(double));
      newParm->bsol = new double[ nphb ];
      memcpy(newParm->bsol, this->bsol, nphb * sizeof(double));
      newParm->hbcut = new double[ nphb ];
      memcpy(newParm->hbcut, this->hbcut, nphb * sizeof(double));
    }
  }  
  // Bond force constant and eq values
  newParm->numbnd = this->numbnd;
  if (this->numbnd > 0) {
    newParm->bond_rk = new double[ numbnd ];
    memcpy(newParm->bond_rk, this->bond_rk, numbnd * sizeof(double));
    newParm->bond_req = new double[ numbnd ];
    memcpy(newParm->bond_req, this->bond_req, numbnd * sizeof(double));
  }
  // Angle force constant and eq values
  newParm->numang = this->numang;
  if (this->numang > 0) {
    newParm->angle_tk = new double[ numang ];
    memcpy(newParm->angle_tk, this->angle_tk, numang * sizeof(double));
    newParm->angle_teq = new double[ numang ];
    memcpy(newParm->angle_teq, this->angle_teq, numang * sizeof(double));
  }
  // Dihedral constant, periodicity, phase
  newParm->numdih = this->numdih;
  if (this->numdih>0) {
    newParm->dihedral_pk = new double[ numdih ];
    memcpy(newParm->dihedral_pk, this->dihedral_pk, numdih * sizeof(double));
    newParm->dihedral_pn = new double[ numdih ];
    memcpy(newParm->dihedral_pn, this->dihedral_pn, numdih * sizeof(double));
    newParm->dihedral_phase = new double[ numdih ];
    memcpy(newParm->dihedral_phase, this->dihedral_phase, numdih * sizeof(double));
  }
  // SCEE and SCNB scale factors
  if (this->scee_scale!=NULL) {
    newParm->scee_scale = new double[ numdih ];
    memcpy(newParm->scee_scale, this->scee_scale, numdih * sizeof(double));
  }
  if (this->scnb_scale!=NULL) {
    newParm->scnb_scale = new double[ numdih ];
    memcpy(newParm->scnb_scale, this->scnb_scale, numdih * sizeof(double));
  }

  // SOLTY is currently unused, just make a straight up copy for now.
  newParm->natyp = this->natyp;
  if (this->natyp > 0) {
    newParm->solty = new double[ natyp ];
    memcpy(newParm->solty, this->solty, natyp * sizeof(double));
  }

  // Fix up IPRES
  newParm->resnums[jres+1] = j;

  // Set up new parm information
  newParm->natom = j;
  newParm->nres = jres+1; 
  newParm->parmFrames = this->parmFrames;
  if (this->molecules>0) 
    newParm->molecules = jmol+1;

  // Set up excluded atoms list
  newParm->SetupExcludedAtoms();

  // Give stripped parm the same pindex as original parm
  newParm->pindex = this->pindex;
  
  // Reallocate memory 
  int *tempresnums = new int[ newParm->nres + 1 ];
  memcpy(tempresnums, newParm->resnums, (newParm->nres + 1) * sizeof(int));
  delete[] newParm->resnums;
  newParm->resnums = tempresnums;
  NAME *tempresnames = new NAME[ newParm->nres ];
  memcpy(tempresnames, newParm->resnames, newParm->nres * sizeof(NAME));
  delete[] newParm->resnames;
  newParm->resnames = tempresnames;
  if (newParm->molecules>0) {
    int *tempAPM = new int[ newParm->molecules ];
    memcpy(tempAPM, newParm->atomsPerMol, newParm->molecules * sizeof(int));
    delete[] newParm->atomsPerMol;
    newParm->atomsPerMol = tempAPM;
  }

  // Set up solvent info if necessary
  if (newParm->firstSolvMol < 0) {
    // No solvent in stripped parmtop
    newParm->solventMolecules=0;
    newParm->hasSolventInfo = false;
  } else {
    // Set up new solvent info based on new resnums and firstSolvMol
    if (newParm->SetSolventInfo()) {
      delete newParm;
      return NULL;
    }
  }
  
  // Copy box information
  for (i=0; i<6; i++)
    newParm->Box[i] = this->Box[i];
  newParm->boxType=this->boxType;


  return newParm;
}

