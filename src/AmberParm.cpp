/* AmberParm.cpp
 * Class that holds parameter information. Can be read in from Amber Topology,
 * PDB, or Mol2 files (implemented in the ReadParmXXX functions). The following
 * parameters of AmberParm must always be set:
 *   The names, resnames, resnums arrays.
 *   The natom, boxType and nres variables.
 * Compiler Defines:
 * - USE_CHARBUFFER: Use CharBuffer to buffer entire file (PDB, Amber Topology) 
 */
#include <cstring>
#include <cstdio> // For sscanf, sprintf
#include <ctime> // for writing time/date to amber parmtop
#include "AmberParm.h" // CpptrajFile.h
#include "PtrajMask.h" // parseMaskString
#include "FortranFormat.h" 
#include "PDBfileRoutines.h"
#include "Mol2FileRoutines.h"
#include "CpptrajStdio.h"
// For searching for bonds by distance (PDB etc)
#include "DistRoutines.h"

#define AMBERPOINTERS 31
#define ELECTOAMBER 18.2223
#define AMBERTOELEC 1/ELECTOAMBER
// =============================================================================

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

// -----------------------------------------------------------------------------
// AmberParm::SetupAtomMask()
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
int AmberParm::SetupIntegerMask(AtomMask &atommaskIn, double *Xin) {
  return SetupAtomMask(atommaskIn, Xin, false);
}

// AmberParm::SetupCharMask()
int AmberParm::SetupCharMask(AtomMask &atommaskIn, double *Xin) {
  return SetupAtomMask(atommaskIn, Xin, true);
}

// -----------------------------------------------------------------------------
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
/** Given a residue number and an atom name, return the atom number. If
  * the given atom name is not in the given residue, return -1.
  */
int AmberParm::FindAtomInResidue(int res, char *atname) {
  if (res < 0 || res >= nres) return -1;
  for (int atnum = resnums[res]; atnum < resnums[res+1]; atnum++) {
    if (strcmp(names[atnum],atname)==0) return atnum;
  }
  return -1;
}

// AmberParm::ResAtomRange()
/// Set the first and last+1 atoms for the given residue
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

// -------------------- ROUTINES FOR ACCESSING INTERNAL DATA -------------------
int AmberParm::NumExcludedAtoms(int atom) {
  if (numex==NULL) return -1;
  if (atom<0 || atom>=natom) return -1;
  return numex[atom];
}

int AmberParm::Natex(int idx) {
  if (excludedAtoms==NULL) return -1;
  return excludedAtoms[idx];
}

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
/// Return bond distance for the two given atoms based on their names.
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
/// Set the atomic charges from the given array.
int AmberParm::SetCharges(double *chargeIn) {
  if (chargeIn==NULL) return 1;
  if (charge==NULL) charge = new double[ natom ];
  memcpy(charge, chargeIn, natom * sizeof(double));
  return 0;
}

// AmberParm::AmberChargeArray()
/// Set the input array with atom charges in Amber format (pre-multiplied by 18.2223)
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
  atom_charge.resize(natom, 18.2223);
  for (int atom = 0; atom < natom; atom++)
    atom_charge[atom] *= charge[atom];
  return 0;
}

// AmberParm::AtomCharge()
/// Return charge on given atom
double AmberParm::AtomCharge(int atomIn) {
  if (charge==NULL) return 0;
  if (atomIn<0 || atomIn >= natom) return 0;
  return charge[atomIn];
}

// AmberParm::AtomsPerMol()
/// Return number of atoms in given molecule.
int AmberParm::AtomsPerMol(int molIn) {
  if (atomsPerMol==NULL) return 0;
  if (molIn < 0 || molIn >= molecules) return 0;
  return atomsPerMol[molIn];
}

// -------------------- ROUTINES PERTAINING TO SURFACE AREA --------------------
// AssignLCPO()
/// Assign parameters for LCPO method. All radii are incremented by 1.4 Ang.
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
  char atype[2];
 
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

// -------------------- ROUTINES PERTAINING TO SOLVENT INFO --------------------
// AmberParm::IsSolventResname()
/// Return true if the residue name corresponds to solvent.
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
      } // END if residue is solvent
        //else mprintf(" not solvent.\n"); // DEBUG
    }
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
    
// --------========= ROUTINES PERTAINING TO READING PARAMETERS =========--------
// AmberParm::OpenParm()
/// Attempt to open file and read in parameters.
int AmberParm::OpenParm(char *filename, bool bondsearch, bool molsearch) {
  CpptrajFile parmfile;
  int err=0;

  if ( parmfile.SetupFile(filename,READ,UNKNOWN_FORMAT, UNKNOWN_TYPE,debug) ) 
    return 1;

  // Copy parm filename to parmName. Separate from File.filename in case of stripped parm
  parmName=new char[ strlen(parmfile.basefilename)+1 ];
  strcpy(parmName,parmfile.basefilename);
  parmfileName=new char[ strlen(filename)+1 ]; 
  strcpy(parmfileName,filename);

  if ( parmfile.OpenFile() ) return 1;

  switch (parmfile.fileFormat) {
    case OLDAMBERPARM: err = ReadParmOldAmber(parmfile); break;
    case AMBERPARM   : err = ReadParmAmber(parmfile);    break;
    case PDBFILE     : err = ReadParmPDB(parmfile)  ;    break;
    case MOL2FILE    : err = ReadParmMol2(&parmfile) ;    break;
    case CHARMMPSF   : err = ReadParmPSF(&parmfile)  ;    break;
    default: 
      rprintf("Unknown parameter file type: %s\n",parmfile.filename);
      err=1;
  }

  parmfile.CloseFile();
  if (err>0) {
    mprinterr("Error reading parm file [%s]\n",filename);
    return 1;
  }

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

/*  // Standardize lengths of atom names and residue names. 4 chars, no
  // leading whitespace. Wrap atom names if they start with a digit, e.g.
  // 1CA becomes CA1. Replace asterisks with ', * is reserved for the mask
  // parser.
  // NOTE: It appears the mask parser is OK with names starting with digits,
  //       so dont worry about that for now.
  for (int atom=0; atom < natom; atom++) { 
    PadWithSpaces(names[atom]);
    TrimName(names[atom]);
    //WrapName(names[atom]);
    ReplaceAsterisk(names[atom]);
  }
  for (int res=0; res < nres; res++) {
    PadWithSpaces(resnames[res]); 
    TrimName(resnames[res]);
    ReplaceAsterisk(names[res]);
  }*/

  // Set up bond information if specified and necessary
  if (bondsearch) {
    if (bonds==NULL && bondsh==NULL && parmCoords!=NULL)
      GetBondsFromCoords();
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

// AmberParm::SetParmFromValues()
/** Used by ReadParmAmber and ReadParmOldAmber to set AmberParm variables
  * from the POINTERS section of the parmtop.
  */
void AmberParm::SetParmFromValues(int *values, bool isOld) {
  // Set some commonly used values
  natom=values[NATOM];
  nres=values[NRES];
  //ifbox=values[IFBOX];
  NbondsWithH=values[NBONH];
  NbondsWithoutH=values[NBONA];
  if (debug>0) {
    if (isOld)
      mprintf("    Old Amber top");
    else
      mprintf("    Amber top");
    mprintf("contains %i atoms, %i residues.\n",natom,nres);
    mprintf("    %i bonds to hydrogen, %i other bonds.\n",NbondsWithH,NbondsWithoutH);
  }
  // Other values
  ntypes = values[NTYPES];
  nnb = values[NNB];
  numbnd = values[NUMBND];
  numang = values[NUMANG];
  numdih = values[NPTRA];
  NanglesWithH=values[NTHETH];
  NanglesWithoutH=values[NTHETA];
  NdihedralsWithH=values[NPHIH];
  NdihedralsWithoutH=values[NPHIA];
  natyp = values[NATYP];
  nphb = values[NPHB];
  // Check that NBONA == MBONA etc. If not print a warning
  if (values[MBONA] != values[NBONA])
    mprintf("\tWarning: [%s] Amber parm has constraint bonds, but they will be ignored.\n");
  if (values[MTHETA] != values[NTHETA])
    mprintf("\tWarning: [%s] Amber parm has constraint angles, but they will be ignored.\n");
  if (values[MPHIA] != values[NPHIA])
    mprintf("\tWarning: [%s] Amber parm has constraint dihedrals, but they will be ignored.\n");
}

// AmberParm::ReadParmOldAmber()
/// Read parameters from an old style (Amber < v7) topology file.
int AmberParm::ReadParmOldAmber(CpptrajFile &parmfile) {
  char *title;
  int values[30], ifbox;

  // TEST: Close and reopen buffered.
  parmfile.CloseFile();
  parmfile.OpenFileBuffered();

  if (debug>0) mprintf("Reading Old-style Amber Topology file %s\n",parmName);
  title = F_load20a4(parmfile);
  if (debug>0) mprintf("\tOld AmberParm Title: %s\n",title);
  delete[] title;
  // Pointers - same as new format except only 30 values, no NEXTRA
  int *tempvalues = (int*) F_loadFormat(parmfile, FINT, 6, 12, 30, debug);
  if (tempvalues==NULL) {
    mprintf("Could not get values from topfile\n");
    return 1;
  }
  memcpy(values, tempvalues, 30 * sizeof(int));
  delete[] tempvalues;
  // Set some commonly used values
  SetParmFromValues(values, true);
  ifbox=values[IFBOX];
  // Load the rest of the parm
  // NOTE: Add error checking!
  names = (NAME*) F_loadFormat(parmfile, FCHAR, 4, 20, natom, debug);
  charge = (double*) F_loadFormat(parmfile, FDOUBLE, 16, 5, natom, debug);
  mass = (double*) F_loadFormat(parmfile, FDOUBLE, 16, 5, natom, debug);
  atype_index = (int*) F_loadFormat(parmfile,FINT, 6, 12, natom, debug);
  numex = (int*) F_loadFormat(parmfile,FINT, 6, 12, natom, debug);
  NB_index = (int*) F_loadFormat(parmfile,FINT, 6, 12, ntypes*ntypes, debug);
  resnames = (NAME*) F_loadFormat(parmfile, FCHAR, 4, 20, nres, debug);
  resnums = (int*) F_loadFormat(parmfile,FINT, 6, 12, nres, debug);
  // Atom #s in resnums are currently shifted +1. Shift back to be consistent
  // with the rest of cpptraj.
  for (int atom=0; atom < nres; atom++)
    resnums[atom] -= 1;
  // Bond, angle, dihedral constants and values
  bond_rk = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NUMBND],debug);
  bond_req = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NUMBND],debug);
  angle_tk = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NUMANG],debug);
  angle_teq = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NUMANG],debug);
  dihedral_pk = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NPTRA],debug);
  dihedral_pn = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NPTRA],debug);
  dihedral_phase = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NPTRA],debug);
  solty = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NATYP],debug);
  // LJ params
  LJ_A = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,ntypes*(ntypes+1)/2,debug);
  LJ_B = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,ntypes*(ntypes+1)/2,debug);
  // Bonds, angles, dihedrals
  bondsh = (int*) F_loadFormat(parmfile,FINT,6,12,values[NBONH]*3,debug);
  bonds = (int*) F_loadFormat(parmfile,FINT,6,12,values[NBONA]*3,debug);
  anglesh = (int*) F_loadFormat(parmfile,FINT,6,12,values[NTHETH]*4,debug); 
  angles = (int*) F_loadFormat(parmfile,FINT,6,12,values[NTHETA]*4,debug); 
  dihedralsh = (int*) F_loadFormat(parmfile,FINT,6,12,values[NPHIH]*5,debug);
  dihedrals = (int*) F_loadFormat(parmfile,FINT,6,12,values[NPHIA]*5,debug);
  // Excluded atoms; shift by -1 so atom #s start from 0
  excludedAtoms = (int*) F_loadFormat(parmfile,FINT,6,12,nnb,debug);
  for (int atom=0; atom < nnb; atom++)
    excludedAtoms[atom] -= 1;
  // LJ 10-12 stuff 
  asol = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NPHB],debug);
  bsol = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NPHB],debug);
  hbcut = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NPHB],debug);
  // Atom types
  types = (NAME*) F_loadFormat(parmfile,FCHAR,4,20,natom,debug);
  // Tree, join, irotat 
  itree = (NAME*) F_loadFormat(parmfile,FCHAR,4,20,natom,debug);
  join_array = (int*) F_loadFormat(parmfile,FINT,6,12,natom,debug);
  irotat = (int*) F_loadFormat(parmfile,FINT,6,12,natom,debug);
  // Solvent/Box info
  if (ifbox > 0) {
    int *solvent_pointer=(int*) F_loadFormat(parmfile,FINT,6,12,3,debug);
    if (solvent_pointer==NULL) {
      mprintf("Error in solvent pointers.\n");
      return 1;
    } else {
      finalSoluteRes=solvent_pointer[0];
      molecules=solvent_pointer[1];
      firstSolvMol=solvent_pointer[2];
      delete[] solvent_pointer;
    }
    atomsPerMol=(int*) F_loadFormat(parmfile,FINT,6,12,molecules,debug);
    if (atomsPerMol==NULL) {mprintf("Error in atoms per molecule.\n"); return 1;}
    // boxFromParm = {OLDBETA, BOX(1), BOX(2), BOX(3)}
    double *boxFromParm=(double*)  F_loadFormat(parmfile,FDOUBLE,16,5,4,debug);
    if (boxFromParm==NULL) {mprintf("Error in box info.\n"); return 1;}
    boxType = SetBoxInfo(boxFromParm,Box,debug);
    delete[] boxFromParm;
    if (debug>0) {
      mprintf("\t%s contains box info: %i mols, first solvent mol is %i\n",
              parmName, molecules, firstSolvMol);
      mprintf("\tBOX: %lf %lf %lf | %lf %lf %lf\n",Box[0],Box[1],Box[2],Box[3],Box[4],Box[5]);
      if (boxType==ORTHO)
        mprintf("\t     Box is orthogonal.\n");
      else if (boxType==NONORTHO)
        mprintf("\t     Box is non-orthogonal.\n");
      else
        mprintf("\t     Box will be determined from first associated trajectory.\n");
    } 
  }
  return 0;
}

// AmberParm::ReadParmAmber() 
/// Read parameters from Amber Topology file
int AmberParm::ReadParmAmber(CpptrajFile &parmfile) {
  int ifbox;
  int *solvent_pointer;
  double *boxFromParm;
  int values[AMBERPOINTERS];
  char *title;
  bool chamber; // true: This topology file is a chamber-created topology file

  // TEST: Close and reopen buffered
  parmfile.CloseFile();
  parmfile.OpenFileBuffered();

  if (debug>0) mprintf("Reading Amber Topology file %s\n",parmName);
  // Title
  // NOTE: getFlagFileString uses 'new' operator.
  title = getFlagFileString(parmfile, "TITLE",debug);
  // If title is NULL, check for CTITLE (chamber parm)
  if (title==NULL) {
    title = getFlagFileString(parmfile,"CTITLE",debug);
    chamber = true;
  } else {
    chamber = false;
  }
  if (debug>0) mprintf("\tAmberParm Title: %s\n",title);
  delete[] title;
  // Pointers
  int *tempvalues=(int*) getFlagFileValues(parmfile,F_POINTERS,AMBERPOINTERS,debug);
  if (tempvalues==NULL) {
    mprintf("Could not get values from topfile\n");
    return 1;
  }
  memcpy(values, tempvalues, AMBERPOINTERS * sizeof(int));
  delete[] tempvalues;
  // Set some commonly used values
  SetParmFromValues(values, false);
  ifbox=values[IFBOX];
  // Atom names
  names=(NAME*) getFlagFileValues(parmfile,F_NAMES,natom,debug);
  if (names==NULL) {mprintf("Error in atom names.\n"); return 1;}
  // Charge; convert to units of electron charge
  charge=(double*) getFlagFileValues(parmfile,F_CHARGE,natom,debug);
  if (charge==NULL) {mprintf("Error in charges.\n"); return 1;}
  for (int atom=0; atom < natom; atom++) charge[atom] *= (AMBERTOELEC);
  // Mass
  mass=(double*) getFlagFileValues(parmfile,F_MASS,natom,debug);
  if (mass==NULL) {mprintf("Error in masses.\n"); return 1;}
  // Atom type index
  atype_index = (int*) getFlagFileValues(parmfile,F_ATYPEIDX,natom,debug);
  if (atype_index==NULL) {mprintf("Error in atom type index.\n"); return 1;}
  // Number of excluded atoms
  numex = (int*) getFlagFileValues(parmfile,F_NUMEX,natom,debug);
  if (numex==NULL) {mprintf("Error in number of excluded atoms.\n"); return 1;}
  // Nonbonded parm index
  NB_index = (int*) getFlagFileValues(parmfile,F_NB_INDEX,ntypes*ntypes,debug);
  if (NB_index==NULL) {mprintf("Error in nonbonded parameter index.\n"); return 1;}
  // Residue names
  resnames=(NAME*) getFlagFileValues(parmfile,F_RESNAMES,nres,debug);
  if (resnames==NULL) {mprintf("Error in residue names.\n"); return 1;}
  // Residue atom #s; shift by -1 so that atom #s start from 0
  resnums=(int*) getFlagFileValues(parmfile,F_RESNUMS,nres,debug);
  if (resnums==NULL) {mprintf("Error in residue numbers.\n"); return 1;}
  for (int res=0; res < nres; res++) resnums[res] -= 1;
  // Bond force constants and equilibrium values
  bond_rk = (double*) getFlagFileValues(parmfile, F_BONDRK, values[NUMBND], debug);
  bond_req = (double*) getFlagFileValues(parmfile, F_BONDREQ, values[NUMBND], debug);
  if (bond_rk==NULL || bond_req==NULL) {mprintf("Error in bond constants.\n"); return 1;}
  // Angle force constants and equilibrium values
  angle_tk = (double*) getFlagFileValues(parmfile, F_ANGLETK, values[NUMANG], debug);
  angle_teq = (double*) getFlagFileValues(parmfile, F_ANGLETEQ, values[NUMANG], debug);
  if (angle_tk==NULL || angle_teq==NULL) {mprintf("Error in angle constants.\n"); return 1;}
  // Dihedral force constants, periodicity, and phase values
  dihedral_pk = (double*) getFlagFileValues(parmfile, F_DIHPK, values[NPTRA], debug);
  dihedral_pn = (double*) getFlagFileValues(parmfile, F_DIHPN, values[NPTRA], debug);
  dihedral_phase = (double*) getFlagFileValues(parmfile, F_DIHPHASE, values[NPTRA], debug);
  if (dihedral_pk==NULL || dihedral_pn==NULL || dihedral_phase==NULL) {
    mprintf("Error in dihedral constants.\n"); return 1;
  }
  // SCEE and SCNB scale factors
  scee_scale = (double*) getFlagFileValues(parmfile, F_SCEE, values[NPTRA], debug);
  scnb_scale = (double*) getFlagFileValues(parmfile, F_SCNB, values[NPTRA], debug);
  // SOLTY: currently unused
  solty = (double*) getFlagFileValues(parmfile,F_SOLTY,values[NATYP],debug);
  // Lennard-Jones A/B coefficient
  LJ_A = (double*) getFlagFileValues(parmfile,F_LJ_A,ntypes*(ntypes+1)/2,debug);
  LJ_B = (double*) getFlagFileValues(parmfile,F_LJ_B,ntypes*(ntypes+1)/2,debug);
  if (LJ_A==NULL || LJ_B==NULL) {mprintf("Error reading LJ parameters.\n"); return 1;}
  // Bond information
  bondsh=(int*) getFlagFileValues(parmfile,F_BONDSH,NbondsWithH*3,debug);
  bonds=(int*) getFlagFileValues(parmfile,F_BONDS,NbondsWithoutH*3,debug);
  if (bondsh==NULL || bonds==NULL) {mprintf("Error in bonds.\n"); return 1;}
  // Angle information
  anglesh = (int*) getFlagFileValues(parmfile,F_ANGLESH, values[NTHETH]*4, debug);
  angles  = (int*) getFlagFileValues(parmfile,F_ANGLES , values[NTHETA]*4, debug);
  if (anglesh==NULL || angles==NULL) {mprintf("Error in angles.\n"); return 1;}
  // Dihedral information
  dihedralsh = (int*) getFlagFileValues(parmfile,F_DIHH, values[NPHIH]*5,  debug);
  dihedrals  = (int*) getFlagFileValues(parmfile,F_DIH , values[NPHIA]*5,  debug);
  if (dihedralsh==NULL || dihedrals==NULL) {mprintf("Error in dihedrals.\n"); return 1;}
  // List of excluded atoms; shift by -1 so atom #s start from 0
  excludedAtoms = (int*) getFlagFileValues(parmfile,F_EXCLUDE,nnb,debug);
  if (excludedAtoms==NULL) {mprintf("Error reading list of excluded atoms.\n"); return 1;}
  for (int atom=0; atom < nnb; atom++) excludedAtoms[atom] -= 1;
  // Hbond LJ 10-12 potential terms and cutoff
  asol  = (double*) getFlagFileValues(parmfile,F_ASOL, values[NPHB],debug);
  bsol  = (double*) getFlagFileValues(parmfile,F_BSOL, values[NPHB],debug);
  hbcut = (double*) getFlagFileValues(parmfile,F_HBCUT,values[NPHB],debug);
  // Amber atom types
  types=(NAME*) getFlagFileValues(parmfile,F_TYPES,natom,debug);
  if (types==NULL) {mprintf("Error in atom types.\n"); return 1;}
  // Tree chain classification and joining info 
  itree = (NAME*) getFlagFileValues(parmfile,F_ITREE,natom,debug);
  join_array = (int*) getFlagFileValues(parmfile,F_JOIN,natom,debug);
  // Last atom that would move if atom i was rotated; unused
  irotat = (int*) getFlagFileValues(parmfile,F_IROTAT,natom,debug);
  // Get solvent info if IFBOX>0
  if (values[IFBOX]>0) {
    solvent_pointer=(int*) getFlagFileValues(parmfile,F_SOLVENT_POINTER,3,debug);
    if (solvent_pointer==NULL) {
      mprintf("Error in solvent pointers.\n");
      return 1;
    } else {
      finalSoluteRes=solvent_pointer[0];
      molecules=solvent_pointer[1];
      firstSolvMol=solvent_pointer[2];
      delete[] solvent_pointer;
    }
    atomsPerMol=(int*) getFlagFileValues(parmfile,F_ATOMSPERMOL,molecules,debug);
    if (atomsPerMol==NULL) {mprintf("Error in atoms per molecule.\n"); return 1;}
    // boxFromParm = {OLDBETA, BOX(1), BOX(2), BOX(3)}
    boxFromParm=(double*) getFlagFileValues(parmfile,F_PARMBOX,4,debug);
    // If no box information present in the parm (such as with Chamber prmtops)
    // set the box info if ifbox = 2, otherwise set to NOBOX; the box info will 
    // eventually be set by angles from the first trajectory associated with 
    // this parm.
    if (boxFromParm==NULL) {
      if (not chamber) mprintf("Warning: Prmtop missing Box information.\n");
      // ifbox 2: truncated octahedron for certain
      if (ifbox == 2) {
        boxType = NONORTHO;
        Box[0] = 0.0; Box[1] = 0.0; Box[2] = 0.0;
        Box[3] = TRUNCOCTBETA;
        Box[4] = TRUNCOCTBETA;
        Box[5] = TRUNCOCTBETA;
      } else
        boxType = NOBOX;
    // Determine box type, set Box angles and lengths from beta (boxFromParm[0])
    } else {
      boxType = SetBoxInfo(boxFromParm,Box,debug);
      delete[] boxFromParm;
    }
    if (debug>0) {
      mprintf("\t%s contains box info: %i mols, first solvent mol is %i\n",
              parmName, molecules, firstSolvMol);
      mprintf("\tBOX: %lf %lf %lf | %lf %lf %lf\n",Box[0],Box[1],Box[2],Box[3],Box[4],Box[5]);
      if (boxType==ORTHO)
        mprintf("\t     Box is orthogonal.\n");
      else if (boxType==NONORTHO)
        mprintf("\t     Box is non-orthogonal.\n");
      else
        mprintf("\t     Box will be determined from first associated trajectory.\n");
    }
  }
  // GB parameters; radius set, radii, and screening parameters
  radius_set = getFlagFileString(parmfile,"RADIUS_SET",debug);
  if (radius_set!=NULL) {
    /*radius_set.assign(title);
    delete[] title;
    // Remove whitespace from GB radius set
    radius_set.erase(std::remove(radius_set.begin(), radius_set.end(), ' '), radius_set.end());*/
    if (debug>0) mprintf("\tRadius Set: %s\n",radius_set);
  }
  gb_radii = (double*) getFlagFileValues(parmfile,F_RADII,natom,debug);
  gb_screen = (double*) getFlagFileValues(parmfile,F_SCREEN,natom,debug);
  if (gb_radii==NULL || gb_screen==NULL) {mprintf("Error reading gb parameters.\n"); return 1;}

  // If parm contains IFCAP or IFPERT info, print a warning since cpptraj
  // currently does not read these in.
  if (values[IFCAP] > 0) 
    mprintf("\tWarning: Parm [%s] contains CAP information, which Cpptraj ignores.\n");
  if (values[IFPERT] > 0)
    mprintf("\tWarning: Parm [%s] contains PERT information, which Cpptraj ignores.\n");

  return 0;
}

// AmberParm::SetAtomsPerMolPDB()
/** Use in ReadParmPDB only, when TER is encountered or end of PDB file
  * update the atomsPerMol array. Take number of atoms in the molecule
  * (calcd as current #atoms - #atoms in previous molecule) as input. 
  * Check if the last residue is solvent; if so, set up solvent information.
  * \return the current number of atoms.
  */
int AmberParm::SetAtomsPerMolPDB(int numAtoms) {
  //mprintf("DEBUG:\tCalling SetAtomsPerMolPDB with %i\n",numAtoms);
  if (numAtoms<1) return 0;
  // Check if the current residue is a solvent molecule
  //mprintf("DEBUG: Checking if %s is solvent.\n",resnames[nres-1]);
  //if (nres>0 && IsSolventResname(resnames[nres-1])) {
  //  if (firstSolvMol==-1) {
  //    firstSolvMol = molecules + 1; // +1 to be consistent w/ Amber top
  //    finalSoluteRes = nres - 1;    // +1 to be consistent w/ Amber top
  //  }
  //}
  int *tempAPM = new int[ molecules+1 ];
  memcpy(tempAPM, atomsPerMol, molecules * sizeof(int));
  delete[] atomsPerMol;
  atomsPerMol = tempAPM;
  atomsPerMol[molecules] = numAtoms;
  molecules++;
  return natom;
}

// AmberParm::ReadParmPDB()
/** Open the PDB file specified by filename and set up topology data.
  * Mask selection requires natom, nres, names, resnames, resnums.
  */
int AmberParm::ReadParmPDB(CpptrajFile &parmfile) {
  char buffer[256];
  int bufferLen;  
  int currResnum;
  int atom;
  int atomInLastMol = 0;
  unsigned int crdidx = 0;

#ifdef USE_CHARBUFFER
  // TEST: Close and reopen buffered.
  parmfile.CloseFile();
  parmfile.OpenFileBuffered();
#endif

  mprintf("    Reading PDB file %s as topology file.\n",parmName);
  currResnum=-1;
  memset(buffer,' ',256);
#ifdef USE_CHARBUFFER
  while ( parmfile.Gets(buffer,256) == 0 ) 
#else
  while ( parmfile.IO->Gets(buffer,256)==0 ) 
#endif
  {
    // If ENDMDL or END is reached stop reading
    if ( strncmp(buffer,"END",3)==0) break;
    // If TER increment number of molecules and continue
    if ( strncmp(buffer,"TER",3)==0) {
      atomInLastMol = SetAtomsPerMolPDB(natom - atomInLastMol);
      continue;
    }
    // Skip all other non-ATOM records
    if (isPDBatomKeyword(buffer)) {
      // Detect and remove trailing newline
      bufferLen = strlen(buffer);
      if (buffer[bufferLen-1] == '\n') buffer[bufferLen-1]='\0';

      // Allocate memory for atom name.
      NAME *tempname = new NAME[ natom+1 ];
      memcpy(tempname, names, natom * sizeof(NAME));
      delete[] names;
      names = tempname;
      // Leading whitespace will automatically be trimmed.
      // Name will be wrapped if it starts with a digit.
      // Asterisks will be replaced with prime char
      pdb_name(buffer, (char*)names[natom]);

      // Allocate memory for coords
      double *tempcoord = new double[ (natom+1)*3 ];
      memcpy(tempcoord, parmCoords, (natom*3)*sizeof(double) );
      delete[] parmCoords;
      parmCoords = tempcoord;
      pdb_xyz(buffer, parmCoords + crdidx);
      crdidx+=3;

      // If this residue number is different than the last, allocate mem for new res
      if (currResnum!=pdb_resnum(buffer)) {
        NAME *temprname = new NAME[ nres+1 ];
        memcpy(temprname, resnames, nres * sizeof(NAME));
        delete[] resnames;
        resnames = temprname;
        // Leading whitespace will automatically be trimmed.
        // Asterisks will be replaced with prime char
        pdb_resname(buffer, (char*)resnames[nres]);
        if (debug>3) mprintf("        PDBRes %i [%s]\n",nres,resnames[nres]);
        int *temprnum = new int[ nres+1 ];
        memcpy(temprnum, resnums, nres * sizeof(int));
        delete[] resnums;
        resnums = temprnum;
        resnums[nres]=natom; 
        currResnum=pdb_resnum(buffer);
        nres++;
        
        // If new residue and HETATM consider it a different molecule as well
        if (strncmp(buffer,"HETATM",6)==0) {
          // If HETATM immediately preceded by a TER card atomsPerMol has
          // just been set, so would be calling with 0. No need to call.
          if ( (natom - atomInLastMol) != 0)
            atomInLastMol = SetAtomsPerMolPDB(natom - atomInLastMol);
        }
  
      // If residue number hasnt changed check for duplicate atom names in res
      // NOTE: At this point nres has been incremented. Want nres-1.
      //       natom is the current atom.
      } else {
        for (atom=resnums[nres-1]; atom < natom; atom++) {
          if ( strcmp(names[natom], names[atom])==0 ) {
            mprintf("      Warning: Duplicate atom name in residue %i [%s]:%i\n",
                    nres,names[natom],natom+1);
          }
        }
      }
      // Clear the buffer
      memset(buffer,' ',256);

      natom++;
    } // END if atom/hetatm keyword
  } // END read in parmfile

  // If a TER card has been read and we are setting up the number of molecules,
  // finish up info on the last molecule read.
  if (molecules>0) {
    SetAtomsPerMolPDB(natom - atomInLastMol);
    // DEBUG
    if (debug>0) {
      //mprintf("\tPDB: firstSolvMol= %i\n",firstSolvMol);
      mprintf("\tPDB: finalSoluteRes= %i\n",finalSoluteRes);
      if (debug>1) {
        mprintf("\tPDB: Atoms Per Molecule:\n");
        for (atom=0; atom < molecules; atom++) {
          mprintf("\t     %8i %8i\n",atom,atomsPerMol[atom]);
        } 
      }
    }
  }

  // No box for PDB - maybe change later to include unit cell info?
  boxType = NOBOX;

  if (debug>0) 
    mprintf("\tPDB contains %i atoms, %i residues, %i molecules.\n",
            natom,nres,molecules);
  // If no atoms, probably issue with PDB file
  if (natom<=0) {
    mprintf("Error: No atoms in PDB file.\n");
    return 1;
  }

  return 0;
}

// AmberParm::ReadParmMol2()
/// Read file as a Tripos Mol2 file.
int AmberParm::ReadParmMol2(CpptrajFile *parmfile) {
  char buffer[MOL2BUFFERSIZE];
  int mol2bonds;
  int resnum, currentResnum;
  unsigned int crdidx = 0;
  char resName[5];

  currentResnum=-1;
  mprintf("    Reading Mol2 file %s as topology file.\n",parmName);
  // Get @<TRIPOS>MOLECULE information
  if (Mol2ScanTo(parmfile, MOLECULE)) return 1;
  //   Scan title
  if ( parmfile->IO->Gets(buffer,MOL2BUFFERSIZE) ) return 1;
  if (debug>0) mprintf("      Mol2 Title: [%s]\n",buffer);
  //   Scan # atoms and bonds
  // num_atoms [num_bonds [num_subst [num_feat [num_sets]]]]
  if ( parmfile->IO->Gets(buffer,MOL2BUFFERSIZE) ) return 1;
  mol2bonds=0;
  sscanf(buffer,"%i %i",&natom, &mol2bonds);
  if (debug>0) {
    mprintf("      Mol2 #atoms: %i\n",natom);
    mprintf("      Mol2 #bonds: %i\n",mol2bonds);
  }

  // Allocate memory for atom names, types, and charges.
  names = new NAME[ natom ];
  types = new NAME[ natom ];
  charge = new double[ natom ];
  // Allocate space for coords
  parmCoords = new double[ natom * 3 ];

  // Get @<TRIPOS>ATOM information
  if (Mol2ScanTo(parmfile, ATOM)) return 1;
  for (int atom=0; atom < natom; atom++) {
    if ( parmfile->IO->Gets(buffer,MOL2BUFFERSIZE) ) return 1;
    // atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]]]
    //sscanf(buffer,"%*i %s %*f %*f %*f %s %i %s %lf", names[atom], types[atom],
    //       &resnum,resName, charge+atom);
    //mprintf("      %i %s %s %i %s %lf\n",atom,names[atom],types[atom],resnum,resName,charge[atom]);
    Mol2AtomName(buffer,names[atom]);
    Mol2AtomType(buffer,types[atom]);
    Mol2XYZ(buffer,parmCoords + crdidx);
    crdidx += 3;
    Mol2ResNumName(buffer,&resnum,resName);
    charge[atom]=Mol2Charge(buffer);
    // Check if residue number has changed - if so record it
    if (resnum != currentResnum) {
      NAME *temprname = new NAME[ nres+1 ];
      memcpy(temprname, resnames, nres * sizeof(NAME));
      delete[] resnames;
      resnames = temprname;
      strcpy(resnames[nres], resName);
      int *temprnum = new int[ nres+1 ];
      memcpy(temprnum, resnums, nres * sizeof(int));
      delete[] resnums;
      resnums = temprnum;
      resnums[nres]=atom; 
      currentResnum = resnum;
      nres++;
    }
  }

  // Get @<TRIPOS>BOND information [optional]
  NbondsWithoutH=0;
  NbondsWithH=0;
  if (Mol2ScanTo(parmfile, BOND)==0) {
    for (int bond=0; bond < mol2bonds; bond++) {
      if ( parmfile->IO->Gets(buffer,MOL2BUFFERSIZE) ) return 1;
      // bond_id origin_atom_id target_atom_id bond_type [status_bits]
      //         resnum         currentResnum
      sscanf(buffer,"%*i %i %i\n",&resnum,&currentResnum);
      // mol2 atom #s start from 1
      AddBond(resnum-1, currentResnum-1,0);
    }
  } else {
    mprintf("      Mol2 file does not contain bond information.\n");
  }

  // No box
  boxType = NOBOX;

  mprintf("    Mol2 contains %i atoms, %i residues,\n", natom,nres);
  mprintf("    %i bonds to H, %i other bonds.\n", NbondsWithH,NbondsWithoutH);

  return 0;
}

// AmberParm::ReadParmPSF()
/** Open the Charmm PSF file specified by filename and set up topology data.
  * Mask selection requires natom, nres, names, resnames, resnums.
  */
int AmberParm::ReadParmPSF(CpptrajFile *parmfile) {
  char buffer[256],tag[256];
  NAME psfname, psfresname;
  int bondatoms[8];
  int currResnum;
  int psfresnum;
  int psfattype;
  int nbond,nlines;

  mprintf("    Reading Charmm PSF file %s as topology file.\n",parmName);
  currResnum=-1;
  memset(buffer,' ',256);
  memset(tag,' ',256);
  tag[0]='\0';

  // Read the first line, should contain PSF...
  if (parmfile->IO->Gets(buffer,256)) return 1;
  // Sanity check
  if (buffer[0]!='P' || buffer[1]!='S' || buffer[2]!='F') {
    mprinterr("Error: ReadParmPSF(): Could not read Charmm PSF file.\n");
    return 1;
  }
  // Advance to <natom> !NATOM
  while (strncmp(tag,"!NATOM",6)!=0) {
    if (parmfile->IO->Gets(buffer,256)) return 1;
    sscanf(buffer,"%i %s",&natom,tag);
  }
  mprintf("\tPSF: !NATOM tag found, natom=%i\n",natom);
  // If no atoms, probably issue with PSF file
  if (natom<=0) {
    mprintf("Error: No atoms in PSF file.\n");
    return 1;
  }

  // Allocate memory for atom name, charge, mass.
  names=(NAME*)    new NAME[ natom ];
  mass=(double*)   new double[ natom ];
  charge=(double*) new double[ natom ];

  // Read the next natom lines
  for (int atom=0; atom < natom; atom++) {
    if (parmfile->IO->Gets(buffer,256) ) {
      mprinterr("Error: ReadParmPSF(): Reading atom %i\n",atom+1);
      return 1;
    }
    // Detect and remove trailing newline
    //bufferLen = strlen(buffer);
    //if (buffer[bufferLen-1] == '\n') buffer[bufferLen-1]='\0';
    // Read line
    // ATOM# SEGID RES# RES ATNAME ATTYPE CHRG MASS (REST OF COLUMNS ARE LIKELY FOR CMAP AND CHEQ)
    sscanf(buffer,"%*8i %*4s %i %4s %4s %4i %14lf %14lf",&psfresnum,psfresname,psfname,
           &psfattype,charge+atom,mass+atom);
    // Ensure name has 4 chars
    PadWithSpaces( psfname );
    strcpy(names[atom],psfname);
    // If this residue number is different than the last, allocate mem for new res
    if (currResnum!=psfresnum) {
        NAME *temprname = new NAME[ nres+1 ];
        memcpy(temprname, resnames, nres * sizeof(NAME));
        delete[] resnames;
        resnames = temprname;
        // Ensure resname has 4 chars
        PadWithSpaces( psfresname );
        strcpy(resnames[nres],psfresname);
        if (debug>3) mprintf("        PSFRes %i [%s]\n",nres,resnames[nres]);
        int *temprnum = new int[ nres+1 ];
        memcpy(temprnum, resnums, nres * sizeof(int));
        delete[] resnums;
        resnums = temprnum;
        resnums[nres]=atom; 
        currResnum=psfresnum;
        nres++;
    }
    // Clear the buffer
    memset(buffer,' ',256);
  } // END loop over atoms 

  // Advance to <nbond> !NBOND
  while (strncmp(tag,"!NBOND",6)!=0) {
    if (parmfile->IO->Gets(buffer,256)) return 1;
    sscanf(buffer,"%i %s",&nbond,tag);
  }
  nlines = nbond / 4;
  if ( (nbond % 4) != 0) nlines++;
  for (int bondline=0; bondline < nlines; bondline++) {
    if (parmfile->IO->Gets(buffer,256) ) {
      mprinterr("Error: ReadParmPSF(): Reading bond line %i\n",bondline+1);
      return 1;
    }
    // Each line has 4 pairs of atom numbers
    int nbondsread = sscanf(buffer,"%i %i %i %i %i %i %i %i",bondatoms,bondatoms+1,
                            bondatoms+2,bondatoms+3, bondatoms+4,bondatoms+5,
                            bondatoms+6,bondatoms+7);
    // NOTE: Charmm atom nums start from 1
    for (int bondidx=0; bondidx < nbondsread; bondidx+=2)
      AddBond(bondatoms[bondidx]-1,bondatoms[bondidx+1]-1,-1);
  }
  //mprintf("DEBUG: Charmm PSF last line after bond read:\n");
  //mprintf("\t[%s]\n",buffer);
  mprintf("\t%i bonds to hydrogen.\n\t%i bonds to non-hydrogen.\n",NbondsWithH,NbondsWithoutH);
    

  //boxType = NOBOX;

  //if (debug>0) 
    mprintf("    PSF contains %i atoms, %i residues, %i molecules.\n",
            natom,nres,molecules);

  return 0;
}

// -----------------------------------------------------------------------------
// AmberParm::AtomInfo()
/// Print parm information for atom.
void AmberParm::AtomInfo(int atom) {
  int res = atomToResidue(atom);
  mprintf("  Atom %i:",atom+1);
  mprintf("[%s]",names[atom]);
  mprintf(" Res %i:",res+1);
  mprintf("[%s]",resnames[res]);
  mprintf(" Mol %i", atomToMolecule(atom)+1);
  if (types!=NULL)
    mprintf(" Type=[%s]",types[atom]);
  if (charge!=NULL)
    mprintf(" Charge=%lf",charge[atom]);
  if (mass!=NULL)
    mprintf(" Mass=%lf",mass[atom]);
  mprintf("\n");
}

// AmberParm::Info()
/// Print information about this parm to buffer.
void AmberParm::ParmInfo() {
  if (parmfileName!=NULL)
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
  }
}

// AmberParm::PrintBondInfo()
/// Print information contained in bonds and bondsh arrays.
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
/// Print information on molecules in PRMTOP
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
/// Print information on residues in PRMTOP
void AmberParm::PrintResidueInfo() {
  mprintf("RESIDUES:\n");
  for (int res = 0; res < nres; res++) 
    mprintf("\tResidue %i, %s, first atom %i\n",res+1,resnames[res],resnums[res]+1);
}

// -----------------------------------------------------------------------------
// AmberParm::atomToResidue()
/// Given an atom number, return corresponding residue number.
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
/// Given an atom number, return corresponding molecule number.
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
/// Given an atom number, return corresponding solvent molecule
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

// -----------------------------------------------------------------------------
// AmberParm::ResetBondInfo()
/// Reset the bonds and bondsh arrays, as well as NBONH and NBONA
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
void AmberParm::GetBondsFromCoords() {
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

  if (parmCoords==NULL) return;
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
        D2 = DIST2_NoImage(parmCoords + idx1, parmCoords + idx2);
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
        D2 = DIST2_NoImage(parmCoords + idx1, parmCoords + idx2);
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

/// Setup excluded atoms list and numex based on bond info
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

// -----------------------------------------------------------------------------
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
  int *newatoms;
  int oldNX, newNX;
  int newatm;
  int Nsequence1;
  int mapatom;
  bool atomIsNegative;
  
  if (atomMap==NULL || oldArray==NULL) return NULL;
  // Actual # entries in oldArray
  oldNX = oldN * Nsequence;
  Nsequence1 = Nsequence - 1;
  // Set initial size of new array to that of old array
  newArray = new int[ oldN * Nsequence ];
  newNX=0;
  // NOTE: Make newatoms static?
  newatoms = new int[ Nsequence ];
  // Go through old array, use atomMap to determine what goes into newParm
  for (int oldi=0; oldi < oldNX; oldi += Nsequence) {
    // Check that atoms 0 to Nsequence exist in newParm. If any of the atoms
    // do not exist in newParm bail.
    newatm = -1;
    for (int sequencei = 0; sequencei < Nsequence1; sequencei++) {
      int arrayAtom = oldArray[oldi+sequencei];
      // For dihedrals the atom # can be negative. Convert to positive
      // for use in the atom map.
      if (arrayAtom < 0) {
        mapatom = -arrayAtom;
        atomIsNegative = true;
      } else {
        mapatom = arrayAtom;
        atomIsNegative = false;
      }
      newatm = atomMap[ mapatom / 3 ];
      if (newatm == -1) break;
      if (atomIsNegative)
        newatoms[sequencei] = -newatm;
      else
        newatoms[sequencei] = newatm;
    }
    // If newatm is -1 here that means it didnt exist in newParm for this
    // sequence. Skip the entire sequence.
    if (newatm==-1) continue;
    // Store the final number of the sequence, which is an index
    newatoms[Nsequence1] = oldArray[oldi+Nsequence1];
    // Place the atoms in newatoms in newArray
    for (int sequencei = 0; sequencei < Nsequence1; sequencei++) 
      newArray[newNX+sequencei] = newatoms[sequencei] * 3;
    // Place the index in newatoms in newArray
    newArray[newNX+Nsequence1] = newatoms[Nsequence1];
    // Increment new counter
    newNX += Nsequence;
  }
  delete[] newatoms;
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
// Adapted from ptraj
/**  The goal of this routine is to create a new AmberParm (newParm)
  *  based on the current AmberParm (this), deleting atoms that are
  *  not in the Selected array.
  */
// NOTE: Make all solvent/box related info dependent on IFBOX only?
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

  // Set stripped parm name based on prefix: <prefix>.<oldparmname>
  // If no prefix given set name as: strip.<oldparmname>
  if (prefix!=NULL) {
    newParm->parmName=new char[ strlen(this->parmName)+strlen(prefix)+2 ];
    sprintf(newParm->parmName,"%s.%s",prefix,this->parmName);
  } else {
    newParm->parmName=new char[ strlen(this->parmName)+7];
    sprintf(newParm->parmName,"strip.%s",this->parmName);
  }

  return newParm;
}

// -----------------------------------------------------------------------------
// AmberParm::WriteAmberParm()
/// Write out information from current AmberParm to an Amber parm file
int AmberParm::WriteAmberParm(char *filename) {
  CpptrajFile outfile;
  CharBuffer buffer;
  int solvent_pointer[3];
  int values[AMBERPOINTERS];
  double parmBox[4];
  // For date and time
  time_t rawtime;
  struct tm *timeinfo;
  int largestRes=0; // For determining nmxrs
  int *tempResnums = NULL;
  double *tempCharge = NULL;

  if (parmName==NULL) return 1;

  if ( outfile.SetupFile(filename, WRITE, AMBERPARM, STANDARD, debug) )
    return 1;

  if (outfile.OpenFile()) return 1;

  // HEADER AND TITLE (4 lines, version, flag, format, title)
  buffer.Allocate( 324 ); // (81 * 4), no space for NULL needed since using NewLine() 
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  // VERSION
  // NOTE: tm_mon ranges from 0
  buffer.Sprintf("%-44s%02i/%02i/%02i  %02i:%02i:%02i                  \n",
                     "%VERSION  VERSION_STAMP = V0001.000  DATE = ",
                     timeinfo->tm_mon+1,timeinfo->tm_mday,timeinfo->tm_year%100,
                     timeinfo->tm_hour,timeinfo->tm_min,timeinfo->tm_sec);
  // TITLE
  buffer.Sprintf("%-80s\n%-80s\n%-80s","%FLAG TITLE","%FORMAT(20a4)","");
  buffer.NewLine();
  //outfile.IO->Printf("%-80s\n",parmName);

  // Shift atom #s in resnums by +1 to be consistent with AMBER
  // Also determine # atoms in largest residue for nmxrs
  tempResnums = new int[ nres ];
  memcpy(tempResnums, resnums, nres * sizeof(int));
  for (int res=0; res < nres; res++) {
    int diff = resnums[res+1] - resnums[res];
    if (diff > largestRes) largestRes = diff;
    tempResnums[res] += 1;
  }

  // POINTERS
  memset(values, 0, AMBERPOINTERS * sizeof(int));
  values[NATOM]=natom;
  values[NTYPES]=ntypes;
  values[NBONH]=NbondsWithH;
  values[MBONA]=NbondsWithoutH;
  values[NTHETH]=NanglesWithH;
  values[MTHETA]=NanglesWithoutH;
  values[NPHIH]=NdihedralsWithH;
  values[MPHIA]=NdihedralsWithoutH;
  values[NNB]=nnb;
  values[NRES]=nres;
  //   NOTE: Assuming NBONA == MBONA etc
  values[NBONA]=NbondsWithoutH;
  values[NTHETA]=NanglesWithoutH;
  values[NPHIA]=NdihedralsWithoutH;
  values[NUMBND]=numbnd;
  values[NUMANG]=numang;
  values[NPTRA]=numdih;
  values[NATYP]=natyp;
  values[NPHB]=nphb;
  values[IFBOX]=AmberIfbox(Box[4]);
  values[NMXRS]=largestRes;
  DataToFortranBuffer(buffer,F_POINTERS, values, NULL, NULL, AMBERPOINTERS);
  // ATOM NAMES
  DataToFortranBuffer(buffer,F_NAMES, NULL, NULL, names, natom);
  // CHARGE - might be null if read from pdb
  if (charge!=NULL) {
    // Convert charges to AMBER charge units
    tempCharge = new double[ natom ];
    memcpy(tempCharge, charge, natom * sizeof(double));
    for (int atom=0; atom<natom; atom++)
      tempCharge[atom] *= (ELECTOAMBER);
    DataToFortranBuffer(buffer,F_CHARGE, NULL, tempCharge, NULL, natom);
    delete[] tempCharge;
  }
  // MASS - might be null if read from pdb
  if (mass!=NULL)  
    DataToFortranBuffer(buffer,F_MASS, NULL, mass, NULL, natom);
  // ATOM_TYPE_INDEX
  if (atype_index!=NULL)
    DataToFortranBuffer(buffer,F_ATYPEIDX, atype_index, NULL, NULL, natom);
  // NUMBER_EXCLUDED_ATOMS
  if (numex!=NULL)
    DataToFortranBuffer(buffer,F_NUMEX, numex, NULL, NULL, natom);
  // NONBONDED_PARM_INDEX
  if (NB_index!=NULL)
    DataToFortranBuffer(buffer,F_NB_INDEX, NB_index, NULL, NULL, ntypes * ntypes); 
  // RESIDUE LABEL - resnames
  DataToFortranBuffer(buffer,F_RESNAMES, NULL, NULL, resnames, nres);
  // RESIDUE POINTER - resnums, IPRES; tempResnums is shifted +1, see above
  DataToFortranBuffer(buffer,F_RESNUMS, tempResnums, NULL, NULL, nres);
  delete[] tempResnums;
  // BOND_FORCE_CONSTANT and EQUIL VALUES
  if (bond_rk!=NULL)
    DataToFortranBuffer(buffer,F_BONDRK, NULL, bond_rk, NULL, numbnd);
  if (bond_req!=NULL)
    DataToFortranBuffer(buffer,F_BONDREQ, NULL, bond_req, NULL, numbnd);
  // ANGLE FORCE CONSTANT AND EQUIL VALUES
  if (angle_tk!=NULL)
    DataToFortranBuffer(buffer,F_ANGLETK, NULL, angle_tk, NULL, numang);
  if (angle_teq!=NULL)
    DataToFortranBuffer(buffer,F_ANGLETEQ, NULL, angle_teq, NULL, numang);
  // DIHEDRAL CONSTANT, PERIODICITY, PHASE
  if (dihedral_pk!=NULL)
    DataToFortranBuffer(buffer,F_DIHPK, NULL, dihedral_pk, NULL, numdih);
  if (dihedral_pn!=NULL)
    DataToFortranBuffer(buffer,F_DIHPN, NULL, dihedral_pn, NULL, numdih);
  if (dihedral_phase!=NULL)
    DataToFortranBuffer(buffer,F_DIHPHASE, NULL, dihedral_phase, NULL, numdih);
  // SCEE and SCNB scaling factors
  if (scee_scale!=NULL)
    DataToFortranBuffer(buffer,F_SCEE, NULL, scee_scale, NULL, numdih);
  if (scnb_scale!=NULL)
    DataToFortranBuffer(buffer,F_SCNB, NULL, scnb_scale, NULL, numdih);
  // SOLTY
  if (solty!=NULL)
    DataToFortranBuffer(buffer,F_SOLTY, NULL, solty, NULL, natyp);
  // Lennard-Jones A/B
  if (LJ_A!=NULL)
    DataToFortranBuffer(buffer,F_LJ_A, NULL, LJ_A, NULL, ntypes*(ntypes+1)/2);
  if (LJ_B!=NULL)
    DataToFortranBuffer(buffer,F_LJ_B, NULL, LJ_B, NULL, ntypes*(ntypes+1)/2); 
  // BONDS INCLUDING HYDROGEN - might be null if read from pdb
  if (bondsh != NULL) 
    DataToFortranBuffer(buffer,F_BONDSH, bondsh, NULL, NULL, NbondsWithH*3);
  // BONDS WITHOUT HYDROGEN - might be null if read from pdb
  if (bonds!=NULL) 
    DataToFortranBuffer(buffer,F_BONDS, bonds, NULL, NULL, NbondsWithoutH*3);
  // ANGLES INCLUDING HYDROGEN
  if (anglesh!=NULL)
    DataToFortranBuffer(buffer,F_ANGLESH, anglesh, NULL, NULL, NanglesWithH*4);
  // ANGLES WITHOUT HYDROGEN
  if (angles!=NULL)
    DataToFortranBuffer(buffer,F_ANGLES, angles, NULL, NULL, NanglesWithoutH*4);
  // DIHEDRALS INCLUDING HYDROGEN
  if (dihedralsh!=NULL)
    DataToFortranBuffer(buffer,F_DIHH, dihedralsh, NULL, NULL, NdihedralsWithH*5);
  // DIHEDRALS WITHOUT H
  if (dihedrals!=NULL)
    DataToFortranBuffer(buffer,F_DIH, dihedrals, NULL, NULL, NdihedralsWithoutH*5);
  // EXCLUDED ATOMS LIST
  // Shift atom #s in excludedAtoms by +1 to be consistent with AMBER
  if (excludedAtoms!=NULL) {
    int *tempexclude = new int[ nnb ];
    memcpy(tempexclude, excludedAtoms, nnb * sizeof(int));
    for (int atom = 0; atom < nnb; atom++) tempexclude[atom] += 1; 
    DataToFortranBuffer(buffer,F_EXCLUDE, tempexclude, NULL, NULL, nnb);
    delete[] tempexclude;
  }
  // HBOND ACOEFF, BCOEFF, HBCUT
  if (asol!=NULL)
    DataToFortranBuffer(buffer,F_ASOL,NULL,asol,NULL,nphb);
  if (bsol!=NULL)
    DataToFortranBuffer(buffer,F_BSOL,NULL,bsol,NULL,nphb);
  if (hbcut!=NULL)
    DataToFortranBuffer(buffer,F_HBCUT,NULL,hbcut,NULL,nphb);
  // AMBER ATOM TYPE - might be null if read from pdb
  if (types!=NULL) 
    DataToFortranBuffer(buffer,F_TYPES, NULL, NULL, types, natom);
  // TREE CHAIN CLASSIFICATION
  if (itree!=NULL)
    DataToFortranBuffer(buffer,F_ITREE, NULL, NULL, itree, natom);
  // JOIN ARRAY
  if (join_array!=NULL)
    DataToFortranBuffer(buffer,F_JOIN, join_array, NULL, NULL, natom);
  // IROTAT
  if (irotat!=NULL)
    DataToFortranBuffer(buffer,F_IROTAT, irotat, NULL, NULL, natom);
  // Write solvent info if IFBOX>0
  if (values[IFBOX]>0) {
    // SOLVENT POINTERS
    if (firstSolvMol!=-1) {
      solvent_pointer[0]=finalSoluteRes;
      solvent_pointer[1]=molecules;
      solvent_pointer[2]=firstSolvMol;
    } else {
      solvent_pointer[0]=nres;
      solvent_pointer[1]=molecules;
      solvent_pointer[2]=molecules+1;
    }
    DataToFortranBuffer(buffer,F_SOLVENT_POINTER, solvent_pointer, NULL, NULL, 3);
    // ATOMS PER MOLECULE
    if (atomsPerMol!=NULL) 
      DataToFortranBuffer(buffer,F_ATOMSPERMOL, atomsPerMol, NULL, NULL, molecules);
    // BOX DIMENSIONS
    parmBox[0] = Box[4]; // beta
    parmBox[1] = Box[0]; // boxX
    parmBox[2] = Box[1]; // boxY
    parmBox[3] = Box[2]; // boxZ
    DataToFortranBuffer(buffer,F_PARMBOX, NULL, parmBox, NULL, 4);
  }
  // GB RADIUS SET
  if (radius_set!=NULL) {
    buffer.IncreaseSize(243); // 3 * 81
    buffer.Sprintf("%-80s\n%-80s\n%-80s","%FLAG RADIUS_SET","%FORMAT(1a80)",radius_set);
    buffer.NewLine();
  }
  // GB RADII
  if (gb_radii!=NULL)
    DataToFortranBuffer(buffer,F_RADII, NULL, gb_radii, NULL, natom);
  // GB SCREENING PARAMETERS
  if (gb_screen!=NULL)
    DataToFortranBuffer(buffer,F_SCREEN, NULL, gb_screen, NULL, natom);

  // Write buffer to file
  outfile.IO->Write(buffer.Buffer(), sizeof(char), buffer.CurrentSize());
  outfile.CloseFile();

  return 0;
}
