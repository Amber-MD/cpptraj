#ifndef INC_ATOMMASK_H
#define INC_ATOMMASK_H
/* 
 * Class: AtomMask
 * AtomMask is used to hold an array of integers that represent atom numbers
 * of atoms selected based on a mask string. The base mask parser is in 
 * AmberParm since it requires access to ipres, atom names, etc. 
 * Although an array of ints becomes larger than a simple character mask once
 * more than 25% of the system is selected, it tends to be faster than the 
 * character array up until about 80% of the system is selected, at which 
 * point the speed is comparable.
 */
#include "AmberParm.h"
class AtomMask {
//    AmberParm *P;     // Topology that selection is based on
//    int N;            // Current position in Selected
    int debug;
  public:
    bool invertMask;  // If true atoms outside the mask will be selected.
    char *maskString; // String specifying atom selection
    int Nselected;    // Number of selected atoms in mask
    int *Selected;    // Int array of selected atom numbers, 1 for each selected atom

    AtomMask();
//    AtomMask(int);
    ~AtomMask();

    AtomMask *Copy();              // Return a copy of this mask
    void AddAtom(int);             // Add given atom to mask
    void PrintMaskAtoms();         // Print atoms in Selected to line
    bool None();                   // Return true if Nselected==0
    void SetMaskString(char*);     // Set the mask string. If NULL, set * (all)
    int SetupMask(AmberParm*,int); // Set up Selected based on maskString and given parm
};
#endif
