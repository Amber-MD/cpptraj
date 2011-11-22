#include <list>
#include <vector>
#include <cstddef> // NULL
#include "Bonds.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
BondInfo::BondInfo() {
  natom = 0;
  Molecule=NULL;
}

// DESTRUCTOR
BondInfo::~BondInfo() {
  if (Molecule!=NULL) delete[] Molecule;
}

// BondInfo::Setup()
/// Set up for the given number of atoms.
int BondInfo::Setup(int natomIn) {
  natom = natomIn;
  Molecule = new bondinfo[ natom ];
  // Initialize molecule info for all atoms
  for (int atom=0; atom < natom; atom++) {
    Molecule[atom].mol      = -1;
    Molecule[atom].maxbonds =  0;
    Molecule[atom].nbonds   =  0;
    Molecule[atom].bond[0]  = -1;
    Molecule[atom].bond[1]  = -1;
    Molecule[atom].bond[2]  = -1;
    Molecule[atom].bond[3]  = -1;
    Molecule[atom].bond[4]  = -1;
    Molecule[atom].bond[5]  = -1;
    Molecule[atom].bond[6]  = -1;
  }
  return 0;
}

// BondInfo::SetValences()
/// Set max valence for atoms based on the given atom names
// NOTE: Just set all to max for now.
void BondInfo::SetValences(NAME *Name) {
  for (int atom = 0; atom < natom; atom++) {
    //Molecule[atom].maxbonds = MaxValence(Name[atom]);
    Molecule[atom].maxbonds = MAXNUMBONDS;
  }
}

// BondInfo::BondAtoms()
/// Create a bond from atom2 to atom1.
int BondInfo::BondAtoms(int atom1, int atom2) {
  // Add atom2 to atom1
  if (Molecule[atom1].nbonds + 1 > Molecule[atom1].maxbonds) {
    mprintf("Warning: BondAtoms: Valences for atom %i maxed (%i)!\n",atom1+1,
            Molecule[atom1].maxbonds);
    return 1;
  }
  int idx = Molecule[atom1].nbonds;
  Molecule[atom1].bond[idx] = atom2;
  Molecule[atom1].nbonds++;
  return 0;
}

// BondInfo::CreateBond()
/// Create bond between both atoms.
int BondInfo::CreateBond(int atom1, int atom2) {
  int err = 0;
  //mprinterr("DBG:\tCreating bond between atom %i and %i\n",atom1,atom2);
  err += BondAtoms(atom1,atom2);
  err += BondAtoms(atom2,atom1);
  return err;
}

// BondInfo::SetBondsFromAmberArray
/// Set up bonding information from Amber-style bond arrays
void BondInfo::SetBondsFromAmberArray(int *bondsIn, int N) {
  int N3 = N * 3;
  if (bondsIn==NULL) return;
  for (int bond = 0; bond < N3; bond += 3) {
    int atom1 = bondsIn[bond  ] / 3;
    int atom2 = bondsIn[bond+1] / 3;
    CreateBond(atom1,atom2);
  }
}

// BondInfo::PrintBonds()
void BondInfo::PrintBonds() {
  for (int atom=0; atom < natom; atom++) {
    mprintf("\t%8i [%8i]:",atom+1,Molecule[atom].mol);
    for (int bond=0; bond < Molecule[atom].nbonds; bond++) {
      mprintf(" %i",Molecule[atom].bond[bond]+1);
    }
    mprintf("\n");
  }
}

// BondInfo::VisitAtom()
void BondInfo::VisitAtom(int atom, int mol) {
  // If this atom has already been visited return
  if (Molecule[atom].mol!=-1) return;
  // Mark this atom as visited
  Molecule[atom].mol=mol;
  // Visit each atom bonded to this atom
  for (int bond = 0; bond < Molecule[atom].nbonds; bond++) 
    VisitAtom(Molecule[atom].bond[bond], mol);
}

// BondInfo::DetermineMolecules() 
int *BondInfo::DetermineMolecules(int *molecules) {
  int mol=0;
  int *atomsPerMol = NULL;
  // First perform recursive search along bonds to determine molecules
  for (int atom = 0; atom < natom; atom++) {
    if (Molecule[atom].mol==-1) {
      //mprintf("\t\tStarting search for molecule %i at atom %i\n",mol,atom+1);
      VisitAtom(atom, mol);
      mol++;
    }
  }
  mprintf("\t%i molecules.\n",mol);

  // Second count how many atoms are in each molecule
  atomsPerMol = new int[mol];
  for (int molecule = 0; molecule < mol; molecule++)
    atomsPerMol[ molecule ] = 0;
  for (int atom=0; atom < natom; atom++) {
    int molecule = Molecule[atom].mol;
    if (molecule>-1)
      atomsPerMol[ molecule ]++;
  }
  // DEBUG
  //for (int molecule = 0; molecule < mol; molecule++) 
  //  mprintf("\tAtomsPerMol %8i: %i\n",molecule,atomsPerMol[ molecule ]);

  *molecules = mol;
  return atomsPerMol;
}

// BondInfo::AtomDistance()
void BondInfo::AtomDistance(int atom, int dist) {
  // If this atom is already too far away return
  if (dist==4) return;
  // Mark distance for this atom 
  Molecule[atom].mol=dist;
  // Visit each atom bonded to this atom
  for (int bond = 0; bond < Molecule[atom].nbonds; bond++) 
    AtomDistance(Molecule[atom].bond[bond], dist+1);
}

// BondInfo::DetermineExcludedAtoms()
/** For each atom, first perform a truncated recursive search to determine
  * which atoms are within 4 bonds away. Search bails once it is farther 
  * than that. Then set up an ordered list of those atoms, which will then
  * be the excluded atoms list. Also set up numex, the number of excluded
  * atoms for each atom. Numex MUST already be allocated
  */
int *BondInfo::DetermineExcludedAtoms(int *numex, int *nnb) {
  std::list<int> excluded_i;
  std::vector<int> excluded;

  if (numex==NULL) {
    mprinterr("Error: BondInfo::DetermineExcludedAtoms: numex is NULL.\n");
    return NULL;
  }
  // Ensure mol for each atom starts at -1
  for (int atomi = 0; atomi < natom; atomi++)
    Molecule[atomi].mol=-1;

  // Determine excluded atoms for each atom
  for (int atomi = 0; atomi < natom; atomi++) {
    AtomDistance(atomi, 0);
    // Now each atom within 4 bonds has mol set to how far away it is. All 
    // other atoms have -1.
    // Determine which atoms with atom# > this atom are closest. 
    for (int atomj = 0; atomj < natom; atomj++) {
      if (atomj > atomi) {
        if (Molecule[atomj].mol!=-1) {
          excluded_i.push_back(atomj);
        }
      }
      // Reset mol for use with next atomi
      Molecule[atomj].mol = -1;
    }
    // Sort
    excluded_i.sort();
    // If no excluded atoms for this atom, insert -1 as a placeholder
    if (excluded_i.empty())
      excluded_i.push_back(-1);
    // Update numex
    numex[atomi] = (int) excluded_i.size();

    // DEBUG
    //mprintf("DBG: %i: [%u]",atomi+1,excluded_i.size());
    //for (std::list<int>::iterator it = excluded_i.begin(); it!=excluded_i.end(); it++)
    //  mprintf(" %i",(*it)+1);
    //mprintf("\n");
    // END DEBUG
    // Append excluded list for i to overall list
    for (std::list<int>::iterator it = excluded_i.begin(); it!=excluded_i.end(); it++)
      excluded.push_back(*it);
    // Clear excluded list for i 
    excluded_i.clear();
  }
  // Convert excluded to int array
  *nnb = (int) excluded.size();
  int *excludedAtoms = new int[ *nnb ];
  int nex = 0;
  for (std::vector<int>::iterator exatom = excluded.begin();
                                  exatom != excluded.end();
                                  exatom++)
  {
    excludedAtoms[nex++] = *exatom;
  }
  return excludedAtoms;
}

// BondInfo::GetListOfBondedAtoms()
/// Get the list of atoms bonded to the given atom
void BondInfo::GetListOfBondedAtoms(int atomIn, int *bondList, int *nbonds) {
  *nbonds = 0;
  if (atomIn < 0 || atomIn > natom) return;
  for (int bond = 0; bond < Molecule[atomIn].nbonds; bond++)
    bondList[bond] = Molecule[atomIn].bond[bond];
  *nbonds = Molecule[atomIn].nbonds;
}

// BondInfo::MaskOfAtomsAroundBond()
/** Given that atom1 and atom2 are bonded, return a mask consisting
  * of the atom #s of all atoms bonded to atom2, excluding atom1
  * and all atoms bonded to atom1.
  */
int *BondInfo::MaskOfAtomsAroundBond(int atom1, int atom2, int *Nselected) {
  int *Selected = NULL;
  int N = 0;
  *Nselected = 0;

  // Ensure atom1 bonded to atom2?
  // First select atom2
  Molecule[atom2].mol = 1;
  // Mark each atom bonded to atom2 (not including atom1) recursively
  for (int bond = 0; bond < Molecule[atom2].nbonds; bond++) {
    if ( Molecule[atom2].bond[bond] != atom1 ) {
      VisitAtom( Molecule[atom2].bond[bond], 1);
    }
  }
  // Count the number of atoms selected
  for (int atom = 0; atom < natom; atom++) 
    if ( Molecule[atom].mol==1 ) N++;
  //mprintf("DEBUG:\tMaskOfAtomsAroundBond(%i,%i) %i atoms selected.\n",atom1,atom2,N);
  if (N==0) return NULL;
  Selected = new int[ N ];
  // Fill selected array, resetting mol as we go
  N = 0;
  for (int atom = 0; atom < natom; atom++) {
    if ( Molecule[atom].mol==1 ) Selected[N++] = atom;
    Molecule[atom].mol=-1;
  }
  *Nselected = N;
  return Selected;
}

/* ========================================================================== */
/* AtomicNumberFromName()
 */
int AtomicNumberFromName(char *Name) {
  char *ptr;
  int element = -1;
  // position ptr at first non-space character in name
  ptr=Name;
  while (*ptr==' ' && *ptr!='\0') ptr++;
  // if NULL something went wrong, abort
  if (*ptr=='\0') {
    return -1;
  }
  switch (ptr[0]) {
    case 'H' : element = 1; break;
    case 'B' :
      if (ptr[1]=='r' || ptr[1]=='R') element = 35;
      else element = 5;
      break;
    case 'C' :
      if (ptr[1]=='l' || ptr[1]=='L') element = 17;
      else element = 6;
      break;
    case 'N' : element = 7; break;
    case 'O' : element = 8; break;
    case 'F' :
      if (ptr[1]=='e' || ptr[1]=='E') element = 26;
      else element = 9;
      break;
    case 'P' : element = 15; break;
    case 'S' : element = 16; break;
    default:
      mprintf("Warning: Could not determine atomic number from name [%s]\n",Name);
  }
  return element;
}
 
/* GetElementFromName()
 * Base on atom name, return 1 character signifying the element.
 * Convert chlorine to X, bromine to Y
 */
char ElementFromName(char *Name) {
  char *ptr;
  char element = 0;
  // position ptr at first non-space character in name
  ptr=Name;
  while (*ptr==' ' && *ptr!='\0') ptr++;
  // if NULL something went wrong, abort
  if (*ptr=='\0') {
    return element;
  }
  element=ptr[0];
  // If no more characters return now
  if (ptr[1]=='\0') return element;
  // If C, check for L or l for chlorine
  if (ptr[0]=='C') {
    if (ptr[1]=='L' || ptr[1]=='l') element='X';
  }
  // If B, check for R or r for bromine
  if (ptr[0]=='B') {
    if (ptr[1]=='R' || ptr[1]=='r') element='Y';
  }
  // DEBUG
  //mprintf("\t\tAtom %s element: [%c]\n",Name,element);
  return element;
}

/* compareElement()
 * Compare a pair of characters a1 and a2 with reference values b1 and b2.
 * If either combination of a1 and a2 match b1 and b2 (i.e. a1==b1 and 
 * a2==b2, or a1==b2 and a2==b1) return true, otherwise return false.
 */
static bool compareElement(char a1, char a2, const char b1, const char b2) {
  if      (a1==b1 && a2==b2) 
    return true;
  else if (a1==b2 && a2==b1) 
    return true;
  return false;
}

/* MaxValence()
 * For the given atom name, determine the element with ElementFromName
 * and return the maximum number of bonds.
 * NOTE: Although Cl and Br (X and Y) can theoretically have 5 bonds,
 * in most cases they should only have 1, so return 1.
*/
int MaxValence(char *Name) {
  char element = ElementFromName(Name);
  int valence = MAXNUMBONDS; // default
  switch (element) {
    case 'H' : valence = 1; break;
    case 'B' : valence = 3; break;
    case 'C' : valence = 4; break;
    case 'N' : valence = 3; break;
    case 'O' : valence = 2; break;
    case 'F' : valence = 1; break;
    case 'P' : valence = 5; break;
    case 'S' : valence = 6; break;
    case 'X' : valence = 1; break;
    case 'Y' : valence = 1; break;
    default:
      mprintf("\tWarning: Could not determine element for %s, using max valence (%i)\n",
              valence);
  }
  return valence;
}

/* IsIon()
 * Check for a '-' or '+' anywhere after the first char of the name;
 * if found assume this is an ion.
 */
static bool IsIon(char *Name) {
  char *ptr = Name + 1;
  while (*ptr!='\0') {
    if (*ptr=='-' || *ptr=='+') return true;
    ptr++;
  }
  return false;
}

/* GetBondedCut()
 * See below. This version converts atom names to 1 char element codes, then
 * determines the cutoff.
 */
double GetBondedCut(char *A1, char *A2) {
  char atom1;
  char atom2;
  // First check for ions, which will never covalently bond
  if (IsIon(A1)) return 0;
  if (IsIon(A2)) return 0;
  // Convert atom names to element names.
  atom1 = ElementFromName(A1);
  atom2 = ElementFromName(A2);
  return GetBondedCut(atom1,atom2);
}

/* GetBondedCut() 
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
double GetBondedCut(char atom1, char atom2) {
  // Default cutoff
  double cut=1.60;
  // Self
  if (atom1==atom2) {
    if      (atom1=='H') cut=0.74;
    else if (atom1=='N') cut=1.45;
    else if (atom1=='C') cut=1.54;
    else if (atom1=='O') cut=1.48;
    else if (atom1=='P') cut=2.21;
    else if (atom1=='S') cut=2.05; // S-S gas-phase value; S=S is 1.49
  }
  // Bonds to H 
  else if ( compareElement(atom1,atom2,'H','C') )
    cut=1.09;
  else if ( compareElement(atom1,atom2,'H','N') )
    cut=1.01;
  else if ( compareElement(atom1,atom2,'H','O') )
    cut=0.96;
  else if ( compareElement(atom1,atom2,'H','P') )
    cut=1.44;
  else if ( compareElement(atom1,atom2,'H','S') )
    cut=1.34;
  // Bonds to C
  else if ( compareElement(atom1,atom2,'C','N') )
    cut=1.47;
  else if ( compareElement(atom1,atom2,'C','O') )
    cut=1.43;
  else if ( compareElement(atom1,atom2,'C','P') )
    cut=1.84;
  else if ( compareElement(atom1,atom2,'C','F') )
    cut=1.35;
  else if ( compareElement(atom1,atom2,'C','X') )
    cut=1.77;
  else if ( compareElement(atom1,atom2,'C','Y') )
    cut=1.94;
  else if ( compareElement(atom1,atom2,'C','S') )
    cut=1.82;
  // Bonds to N
  else if ( compareElement(atom1,atom2,'N','O') )
    cut=1.40;
  else if ( compareElement(atom1,atom2,'N','S') )
    cut=1.68; // Postma & Vos, Acta Cryst. (1973) B29, 915
  else if ( compareElement(atom1,atom2,'N','F') )
    cut=1.36;
  else if ( compareElement(atom1,atom2,'N','X') )
    cut=1.75;
  else if ( compareElement(atom1,atom2,'N','P') )
    cut=1.71; // Avg over all nX-pX from gaff.dat
  // Bonds to P
  else if ( compareElement(atom1,atom2,'P','O') )
    cut=1.63;
  else if ( compareElement(atom1,atom2,'P','S') )
    cut=1.86;
  else if ( compareElement(atom1,atom2,'P','F') )
    cut=1.54;
  else if ( compareElement(atom1,atom2,'P','X') )
    cut=2.03; 
  // Bonds to O
  else if ( compareElement(atom1,atom2,'O','S') )
    cut=1.48;
  else if ( compareElement(atom1,atom2,'O','F') )
    cut=1.42;
  // Bonds to S
  else if ( compareElement(atom1,atom2,'S','F') )
    cut=1.56;
  else if ( compareElement(atom1,atom2,'S','X') )
    cut=2.07;
  // No cutoff, use default
  else {
    mprintf("Warning: GetBondedCut: Cut not found for %c - %c\n",atom1,atom2);
    mprintf("                       Using default cutoff of %lf\n",cut);
  }

  // Padding value
  cut+=0.1;
  return cut;
}

