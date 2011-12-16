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
/// Used to determine how many unique non-covalently bonded molecules
/// exist within a set of atoms.
class BondInfo {
    // Mol maxbonds #bonds bond1 bond2 bond3 bond4 bond5 bond6 bond7 
    // 0   1        2      3     4     5     6     7     8     9
    struct bondinfo {
      int mol;
      int maxbonds;
      std::vector<int> bond;
    }; 
    std::vector<bondinfo> Molecule;
    int natom;

    int BondAtoms(int,int);
    void VisitAtom(int, int);
    void AtomDistance(int, int);
  public:
    BondInfo();
    ~BondInfo();
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
char GetElementFromName(char *);
double GetBondedCut(NAME, NAME);
double GetBondedCut(char , char );
double GetBondedCut(AtomicElementType, AtomicElementType);
#endif
