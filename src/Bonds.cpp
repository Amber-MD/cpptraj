#include "Bonds.h"
#include "CpptrajStdio.h"
/* GetElementFromName()
 * Base on atom name, return 1 character signifying the element.
 * Convert chlorine to X, bromine to Y
 */
static char ElementFromName(char *Name) {
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
double GetBondedCut(char *A1, char *A2) {
  char atom1;
  char atom2;
  // Default cutoff
  double cut=1.60;
  // Convert atom names to element names.
  atom1 = ElementFromName(A1);
  atom2 = ElementFromName(A2);
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
    mprintf("Warning: GetBondedCut: Cut not found for %s - %s\n",atom1,atom2);
    mprintf("                       Using default cutoff of %lf\n",cut);
  }

  // Padding value
  cut+=0.1;
  return cut;
}

