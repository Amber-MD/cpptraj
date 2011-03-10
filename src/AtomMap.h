#ifndef INC_ATOMMAP_H
#define INC_ATOMMAP_H
// AtomMap
#include "Action.h"

#define ATOMIDLENGTH 26
#define UNIQUELENGTH 131

class atommap {
    struct mapatom {
      int bond[4];               // Holds indices of other bonded atoms
      int nbond;                 // Number of bonds
      int visited;               // 1=Has been visited by traceBonds
      char atomID[ATOMIDLENGTH]; // ID created from this atom name, then bonded atom names 
      char unique[UNIQUELENGTH]; // ID created from this atomID, then bonded atomIDs
      int isUnique;              // 1 if no other unique ID matches this atom
    };
  public:
    mapatom *M;     // Array of atoms, contains bond info etc
    int natom;      // Number of atoms in the map
    char **names;   // Name of each atom
    Frame *F;       // Hold atom coords
    AmberParm *P;   // Hold corresponding parm
    int debug;

    atommap();
    ~atommap();
    int setup();                   // Set up atom map, get atom elmt names, init map
    double getCut(char *, char *); // Determine bond length cutoff between two atoms
    int calcDist();                // Determine bonds between atoms based on distance
    void determineAtomID();        // Give each atom an id based on what atoms are bonded to it
    void printBonds();             // Print bond info and uniqueness for all atoms
    int fourAtoms(int *);
};

class AtomMap : public Action {
    atommap RefMap;
    atommap TargetMap;
  public:
    AtomMap() {}
    ~AtomMap(){}
    int init();
};
#endif
