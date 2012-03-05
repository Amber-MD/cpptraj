#include <list>
#include <vector>
#include <cstddef> // NULL
#include "Bonds.h"
#include "CpptrajStdio.h"

// ========== FUNCTIONS ========================================================
/// Atom names corresponding to AtomicElementType.
// 2 chars + NULL.
const char AtomicElementName[NUM_DEFINED_ELEMENTS][3] = { "??",
  "H",  "B",  "C",  "N", "O",  "F",  "P",  "S", "Cl", "Br", "Fe", "Ca",
  "I", "Mg", "Cu", "Li", "K", "Rb", "Cs", "Zn", "Na"
};

/// Unique 1 char atom names corresponding to AtomicElementType
// For use when atom mapping
const char AtomicElementChar[NUM_DEFINED_ELEMENTS] = { 0,
  'H',  'B',  'C',  'N',  'O', 'F',  'P',  'S',  'X',  'Y',  'f',  'c',
  'I',  'M',  'U',  'L',  'K', 'R',  'E',  'Z',  'n' 
};

// ElementFromName()
AtomicElementType ElementFromName(NAME Name) {
  char *ptr;
  AtomicElementType element = UNKNOWN_ELEMENT;
  // position ptr at first non-space character in name
  ptr=Name;
  while (*ptr==' ' && *ptr!='\0') ptr++;
  // if NULL something went wrong, abort
  if (*ptr=='\0') return element;
  
  switch (ptr[0]) {
    case 'H' : element = HYDROGEN; break;
    case 'C' :
      if (ptr[1]=='l' || ptr[1]=='L') element = CHLORINE;
      else if (ptr[1]=='0') element = CALCIUM;
      else if (ptr[1]=='s') element = CESIUM;
      else if (ptr[1]=='U') element = COPPER;
      else element = CARBON;
      break;
    case 'N' : 
      if (ptr[1]=='a') element = SODIUM;
      else element = NITROGEN; 
      break;
    case 'O' : element = OXYGEN; break;
    case 'F' :
      if (ptr[1]=='e' || ptr[1]=='E') element = IRON;
      else element = FLUORINE;
      break;
    case 'B' :
      if (ptr[1]=='r' || ptr[1]=='R') element = BROMINE;
      else element = BORON;
      break;
    case 'I' : 
      if (ptr[1]=='M') element = CHLORINE;    // IM, Cl- in old FF94
      else if (ptr[1]=='P') element = SODIUM; // IP, Na+ in old FF94
      else element = IODINE; 
      break;
    case 'P' : element = PHOSPHORUS; break;
    case 'S' : element = SULFUR; break;
    case 'M' :
      if (ptr[1]=='G') element = MAGNESIUM; 
      break;
    case 'Z' :
      if (ptr[1]=='n') element = ZINC;
      break;
    case 'L' :
      if (ptr[1]=='i') element = LITHIUM;
      break;
    case 'K' : element = POTASSIUM; break;
    case 'R' :
      if (ptr[1]=='b') element = RUBIDIUM;
      break; 
    default:
      mprintf("Warning: Could not determine atomic number from name [%s]\n",Name);
  }
  return element;
}

// MaxValence()
/** For the given atom name, determine the element with ElementFromName
  * and return the maximum number of bonds.
  * NOTE: Although Cl and Br can theoretically have 5 bonds,
  * in most cases they should have atom most 1, so return 1.
  */
static int MaxValence(NAME Name) {
  AtomicElementType element = ElementFromName(Name);
  int valence = MAXNUMBONDS; // default
  switch (element) {
    case HYDROGEN   : valence = 1; break;
    case BORON      : valence = 3; break;
    case CARBON     : valence = 4; break;
    case NITROGEN   : valence = 4; break; // Allow 4 bonds for N+
    case OXYGEN     : valence = 2; break;
    case FLUORINE   : valence = 1; break;
    case PHOSPHORUS : valence = 5; break;
    case SULFUR     : valence = 6; break;
    case CHLORINE   : valence = 1; break;
    case BROMINE    : valence = 1; break;
    default:
      mprintf("\tWarning: Valence not found for for [%s], using max valence (%i)\n",
              Name, valence);
  }
  return valence;
}

// compareElement()
/** Compare a pair of elements a1 and a2 with target values b1 and b2.
  * If either combination of a1 and a2 match b1 and b2 (i.e. a1==b1 and 
  * a2==b2, or a1==b2 and a2==b1) return true, otherwise return false.
  */
static bool compareElement(AtomicElementType a1, AtomicElementType a2, 
                           AtomicElementType b1, AtomicElementType b2) 
{
  if      (a1==b1 && a2==b2) 
    return true;
  else if (a1==b2 && a2==b1) 
    return true;
  return false;
}


// IsIon()
/** Check for a '-' or '+' anywhere after the first char of the name;
  * if found assume this is an ion.
  */
static bool IsIon(NAME Name) {
  char *ptr = Name + 1;
  while (*ptr!='\0') {
    if (*ptr=='-' || *ptr=='+') return true;
    ptr++;
  }
  return false;
}

// GetBondedCut()
/** See below. This version converts atom names to elements, then
  * determines the cutoff.
  */
double GetBondedCut(NAME A1, NAME A2) {
  AtomicElementType atom1;
  AtomicElementType atom2;
  // First check for ions, which will never covalently bond
  if (IsIon(A1)) return 0;
  if (IsIon(A2)) return 0;
  // Convert atom names to element names.
  atom1 = ElementFromName(A1);
  atom2 = ElementFromName(A2);
  return GetBondedCut(atom1,atom2);
}

// ---------- AtomMap related routines; eventually obsolete --------------------
// ConvertNameToChar()
/** Base on atom name, return 1 character signifying the element.
  */
char ConvertNameToChar(NAME Name) {
  AtomicElementType element;

  element = ElementFromName(Name);
  if (element == UNKNOWN_ELEMENT) return 0;

  //mprintf("\t\tAtom %s element: [%c]\n",Name,AtomicElementChar[element]);

  // Return corresponding char
  return AtomicElementChar[ element ];
}
// -----------------------------------------------------------------------------

// GetBondedCut() 
/** Return a cutoff based on optimal covalent bond distance based on the 
  * identities of atom1 and atom2. When multiple hybridizations are possible
  * the longest possible bond length is used.
  * Unless otherwise noted values taken from:
  * - Huheey, pps. A-21 to A-34; T.L. Cottrell, "The Strengths of Chemical Bonds," 
  *       2nd ed., Butterworths, London, 1958; 
  * - B. deB. Darwent, "National Standard Reference Data Series," National Bureau of Standards, 
  *       No. 31, Washington, DC, 1970; S.W. Benson, J. Chem. Educ., 42, 502 (1965).
  * Can be found on the web at:
  * - http://www.wiredchemist.com/chemistry/data/bond_energies_lengths.html
  */
// NOTE: Store cutoff^2 instead??
double GetBondedCut(AtomicElementType atom1, AtomicElementType atom2) {
  // Default cutoff
  double cut=1.60;
  // Self
  if (atom1==atom2) {
    if      (atom1==HYDROGEN   ) cut=0.74;
    else if (atom1==NITROGEN   ) cut=1.45;
    else if (atom1==CARBON     ) cut=1.54;
    else if (atom1==OXYGEN     ) cut=1.48;
    else if (atom1==PHOSPHORUS ) cut=2.21;
    else if (atom1==SULFUR     ) cut=2.05; // S-S gas-phase value; S=S is 1.49
  }
  // Bonds to H 
  else if ( compareElement(atom1,atom2,HYDROGEN,CARBON) )
    cut=1.09;
  else if ( compareElement(atom1,atom2,HYDROGEN,NITROGEN) )
    cut=1.01;
  else if ( compareElement(atom1,atom2,HYDROGEN,OXYGEN) )
    cut=0.96;
  else if ( compareElement(atom1,atom2,HYDROGEN,PHOSPHORUS) )
    cut=1.44;
  else if ( compareElement(atom1,atom2,HYDROGEN,SULFUR) )
    cut=1.34;
  // Bonds to C
  else if ( compareElement(atom1,atom2,CARBON,NITROGEN) )
    cut=1.47;
  else if ( compareElement(atom1,atom2,CARBON,OXYGEN) )
    cut=1.43;
  else if ( compareElement(atom1,atom2,CARBON,PHOSPHORUS) )
    cut=1.84;
  else if ( compareElement(atom1,atom2,CARBON,FLUORINE) )
    cut=1.35;
  else if ( compareElement(atom1,atom2,CARBON,CHLORINE) )
    cut=1.77;
  else if ( compareElement(atom1,atom2,CARBON,BROMINE) )
    cut=1.94;
  else if ( compareElement(atom1,atom2,CARBON,SULFUR) )
    cut=1.82;
  // Bonds to N
  else if ( compareElement(atom1,atom2,NITROGEN,OXYGEN) )
    cut=1.40;
  else if ( compareElement(atom1,atom2,NITROGEN,SULFUR) )
    cut=1.68; // Postma & Vos, Acta Cryst. (1973) B29, 915
  else if ( compareElement(atom1,atom2,NITROGEN,FLUORINE) )
    cut=1.36;
  else if ( compareElement(atom1,atom2,NITROGEN,CHLORINE) )
    cut=1.75;
  else if ( compareElement(atom1,atom2,NITROGEN,PHOSPHORUS) )
    cut=1.71; // Avg over all nX-pX from gaff.dat
  // Bonds to P
  else if ( compareElement(atom1,atom2,PHOSPHORUS,OXYGEN) )
    cut=1.63;
  else if ( compareElement(atom1,atom2,PHOSPHORUS,SULFUR) )
    cut=1.86;
  else if ( compareElement(atom1,atom2,PHOSPHORUS,FLUORINE) )
    cut=1.54;
  else if ( compareElement(atom1,atom2,PHOSPHORUS,CHLORINE) )
    cut=2.03; 
  // Bonds to O
  else if ( compareElement(atom1,atom2,OXYGEN,SULFUR) )
    cut=1.48;
  else if ( compareElement(atom1,atom2,OXYGEN,FLUORINE) )
    cut=1.42;
  // Bonds to S
  else if ( compareElement(atom1,atom2,SULFUR,FLUORINE) )
    cut=1.56;
  else if ( compareElement(atom1,atom2,SULFUR,CHLORINE) )
    cut=2.07;
  // No cutoff, use default
  else {
    mprintf("Warning: GetBondedCut: Cut not found for %s - %s\n",
            AtomicElementName[atom1],AtomicElementName[atom2]);
    mprintf("                       Using default cutoff of %lf\n",cut);
  }

  // Padding value
  cut+=0.1;
  return cut;
}

// ========== BONDINFO CLASS ===================================================
// CONSTRUCTOR
BondInfo::BondInfo() {
  natom = 0;
}

// BondInfo::Setup()
/** Set up for the given number of atoms. */
int BondInfo::Setup(int natomIn) {
  bondinfo def;

  natom = natomIn;
  Molecule.clear();
  def.mol      = -1;
  def.maxbonds =  0;

  Molecule.resize( natom, def );
  
  return 0;
}

// BondInfo::Reset()
/** Clear all bond information. */
void BondInfo::Reset() {
  Molecule.clear();
}

// BondInfo::HasBeenSetup()
/** Return true if this already contains bond information.
  */
bool BondInfo::HasBeenSetup() {
  if (!Molecule.empty()) return true;
  return false;
}

// BondInfo::SetValences()
/** Set max valence for atoms based on the given atom names.
  */
// NOTE: Just set all to max for now.
void BondInfo::SetValences(NAME *Name) {
  for (int atom = 0; atom < natom; atom++) {
    Molecule[atom].maxbonds = MaxValence(Name[atom]);
    //Molecule[atom].maxbonds = MAXNUMBONDS;
  }
}

// BondInfo::BondAtoms()
/** Create a bond from atom2 to atom1. */
int BondInfo::BondAtoms(int atom1, int atom2) {
  // Add atom2 to atom1
  Molecule[atom1].bond.push_back(atom2);
  if ((int)Molecule[atom1].bond.size() > Molecule[atom1].maxbonds) {
    // NOTE: Supress this warning for now. Can be a false positive if the atom
    // in question is e.g. a TIP3P hydrogen (which gets 2 bonds).
    //mprintf("Warning: BondAtoms: Valences for atom %i maxed (%i)!\n",atom1+1,
    //        Molecule[atom1].maxbonds);
    return 1;
  }

  return 0;
}

// BondInfo::CreateBond()
/** Create bond between both atoms. */
int BondInfo::CreateBond(int atom1, int atom2) {
  int err = 0;
  //mprinterr("DBG:\tCreating bond between atom %i and %i\n",atom1,atom2);
  err += BondAtoms(atom1,atom2);
  err += BondAtoms(atom2,atom1);
  return err;
}

// BondInfo::SetBondsFromAmberArray
/** Set up bonding information from Amber-style bond arrays */
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
    for (std::vector<int>::iterator bondedatom = Molecule[atom].bond.begin();
                                    bondedatom != Molecule[atom].bond.end();
                                    bondedatom++)
    {
      mprintf(" %i",(*bondedatom)+1);
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
  for (std::vector<int>::iterator bondedatom = Molecule[atom].bond.begin();
                                    bondedatom != Molecule[atom].bond.end();
                                    bondedatom++)
    VisitAtom(*bondedatom, mol);
}

// BondInfo::DetermineMolecules() 
int *BondInfo::DetermineMolecules(int *molecules) {
  int mol=0;
  int *atomsPerMol = NULL;
  // Ensure mol for each atom starts at -1
  for (int atomi = 0; atomi < natom; atomi++)
    Molecule[atomi].mol=-1;
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
  for (std::vector<int>::iterator bondedatom = Molecule[atom].bond.begin();
                                    bondedatom != Molecule[atom].bond.end();
                                    bondedatom++)
    AtomDistance(*bondedatom, dist+1);
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
void BondInfo::GetListOfBondedAtoms(int atomIn, std::vector<int> &bondList) {
  if (atomIn < 0 || atomIn > natom) return;
  bondList = Molecule[atomIn].bond;
}

// BondInfo::MaskOfAtomsAroundBond()
/** Given that atom1 and atom2 are bonded, set up a mask consisting
  * of T for all atoms bonded to atom2 (excluding atom1 and all atoms 
  * bonded to atom1), F for everything else.
  * \return 0 on success, 1 on failure.
  */
int BondInfo::MaskOfAtomsAroundBond(int atom1, int atom2, std::vector<char> &Selected) {
  Selected.clear();
  //int N = 0;

  // Ensure mol for each atom starts at -1
  for (int atomi = 0; atomi < natom; atomi++)
    Molecule[atomi].mol=-1;

  // Ensure atom1 bonded to atom2?
  // First select atom2
  Molecule[atom2].mol = 1;
  // Mark each atom bonded to atom2 (not including atom1) recursively
  for (std::vector<int>::iterator bondedatom = Molecule[atom2].bond.begin();
                                  bondedatom != Molecule[atom2].bond.end();
                                  bondedatom++) 
  {
    if ( *bondedatom != atom1 ) {
      VisitAtom( *bondedatom, 1);
    }
  }
  // Count the number of atoms selected
  //for (int atom = 0; atom < natom; atom++) 
  //  if ( Molecule[atom].mol==1 ) N++;
  //mprintf("DEBUG:\tMaskOfAtomsAroundBond(%i,%i) %i atoms selected.\n",atom1,atom2,N);
  //if (N==0) return 1;
  //Selected.reserve( N );
  Selected.reserve( natom );
  // Fill selected array, resetting mol as we go
  //N = 0;
  for (int atom = 0; atom < natom; atom++) {
    if ( Molecule[atom].mol==1 ) 
      //Selected.push_back( atom );
      Selected.push_back( 'T' );
    else
      Selected.push_back( 'F' );
    Molecule[atom].mol=-1;
  }
  return 0;
}

