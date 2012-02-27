#ifndef INC_BONDS_H
#define INC_BONDS_H
#include <vector>
#include "Name.h"
/*! \file Bonds.h
    \brief Classes and routines used for keeping track of and/or defining covalent bonds.
 */
// ---------- DEFINES
#define MAXNUMBONDS 7
// Class: BondInfo
/// Hold information on what atoms each atom is bonded to.
class BondInfo {
    /// Hold bonding information for an atom
    struct bondinfo {
      int mol;               ///< Molecule number of this atom
      int maxbonds;          ///< Maximum number of bonds this atom should have
      std::vector<int> bond; ///< Array of atom #s that this atom is bonded to
    }; 
    /// Hold bonding information for each atom in the system
    std::vector<bondinfo> Molecule;
    /// Total number of atoms.
    int natom;

    int BondAtoms(int,int);
    void VisitAtom(int, int);
    void AtomDistance(int, int);
  public:
    BondInfo();
    int Setup(int);
    void Reset();
    bool HasBeenSetup();
    void SetValences(NAME*);
    int CreateBond(int,int);
    void SetBondsFromAmberArray(int *, int);
    void PrintBonds();
    int *DetermineMolecules(int*);
    int *DetermineExcludedAtoms(int *, int *);
    void GetListOfBondedAtoms(int, std::vector<int>&);
    int MaskOfAtomsAroundBond(int, int, std::vector<char>&);
};
// Other functions
/// Defined type for elements. Based on entries in $AMBERHOME/dat/parm94.dat and parm10.dat
enum AtomicElementType { UNKNOWN_ELEMENT,  
  HYDROGEN,   BORON,     CARBON,   NITROGEN, OXYGEN,    FLUORINE, 
  PHOSPHORUS, SULFUR,    CHLORINE, BROMINE,  IRON,      CALCIUM,  
  IODINE,     MAGNESIUM, COPPER,   LITHIUM,  POTASSIUM, RUBIDIUM, 
  CESIUM,     ZINC,      SODIUM 
};
#define NUM_DEFINED_ELEMENTS 22
AtomicElementType ElementFromName(char *);
char ConvertNameToChar(NAME); // AtomMap only
double GetBondedCut(NAME, NAME);
double GetBondedCut(AtomicElementType, AtomicElementType);
#endif
