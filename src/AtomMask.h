#ifndef INC_ATOMMASK_H
#define INC_ATOMMASK_H
/// Class: AtomMask
/// AtomMask is used to hold an array of integers that represent atom numbers
/// of atoms selected based on a mask string. This class is actually an 
/// interface to the ptraj mask parser written by Viktor Hornak (ptrajmask.cpp).
/// Takes as input a string and an AmberParm class since the basic parser 
/// requires access to ipres, atom names, etc.
/// First the mask string is set via SetMaskString. Then the actual mask can
/// be set up in two ways.
/// 1) Integer mask
///    Although an array of ints becomes larger than a simple character mask 
///    once more than 25% of the system is selected, it tends to be faster 
///    than the character array up until about 80% of the system is selected, 
///    at which point the speed is comparable. This is the default way to use
///    AtomMask and is how most of the routines in the Frame class have been
///    written to use AtomMask.
/// 2) Character mask
///    This is the original way to use the atom mask, useful e.g. when 
///    you need to know what atoms are not selected as well as what atoms
///    are selected. Unlike the integer mask, the character mask is not
///    directly accessible by outside routines, only by the AtomInCharMask
///    routine.
/// *** NOTE ***
/// Currently no checking is done to prevent the mask from being set up as 
/// both an integer and a char mask!!
#include "AmberParm.h"
class AtomMask {
    char *CharMask;   // Char array of atoms, T if selected, F if not.
    int Nchar;        // Number of chars in CharMask array.
  public:
    bool invertMask;  // If true atoms outside the mask will be selected.
    char *maskString; // String specifying atom selection
    int Nselected;    // Number of selected atoms in mask
    int *Selected;    // Int array of selected atom numbers, 1 for each selected atom

    AtomMask();
    ~AtomMask();

    void Reset();                  // Reset atom mask
    AtomMask & operator=(const AtomMask&);
    AtomMask *CopyMask();          // Return a copy of this mask
    void AddAtom(int);             // Add given atom to mask
    void AddAtoms(int *, int);     // Add a list of atoms to mask
    void AddAtomRange(int,int);    // Add minAtom <= atom < maxAtom to mask
    void PrintMaskAtoms();         // Print atoms in Selected to line
    bool None();                   // Return true if Nselected==0
    void SetMaskString(char*);     // Set the mask string. If NULL, set * (all)
    int SetupMask(AmberParm*,int); // Set up Selected based on maskString and given parm
    int SetupCharMask(AmberParm*,int); // Set up CharMask based on maskString and given parm
    int SetupCharMask(AmberParm*,double*,int); // Set up CharMask potentially using coords 
    bool AtomInCharMask(int atom); // True if given atom is T in CharMask
};
#endif
