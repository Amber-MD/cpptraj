#ifndef INC_ATOMMAP_H
#define INC_ATOMMAP_H
// AtomMap
#include "Action.h"

#define ATOMIDLENGTH 26
#define UNIQUELENGTH 131
#define MAXBONDS 6

class atommap {
    struct mapatom {
      int bond[MAXBONDS];        // Holds indices of other bonded atoms
      int nbond;                 // Number of bonds
      bool complete;             // true: This atom an all bonded atoms have been mapped
      char atomID[ATOMIDLENGTH]; // ID created from this atom name, then bonded atom names 
      char unique[UNIQUELENGTH]; // ID created from this atomID, then bonded atomIDs
      int isUnique;              // 1 if no other unique ID matches this atom
    };
    int debug;
  public:
    mapatom *M;     // Array of atoms, contains bond info etc
    int natom;      // Number of atoms in the map
    char **names;   // Name of each atom
    Frame *F;       // Hold atom coords
    AmberParm *P;   // Hold corresponding parm

    atommap();
    ~atommap();
    void SetDebug(int);
    int setup();                   // Set up atom map, get atom elmt names, init map
    double getCut(char *, char *); // Determine bond length cutoff between two atoms
    int calcDist();                // Determine bonds between atoms based on distance
    void determineAtomID();        // Give each atom an id based on what atoms are bonded to it
    void markComplete();           // Mark unique atoms bonded to all unique atoms as complete
    bool BondIsRepeated(int,int);  // True if bonded atom has same name as another bonded atom
    // DEBUG
    //void WriteMol2(char *);
};

class AtomMap : public Action {
    atommap RefMap;
    atommap TargetMap;
    int mapSingleBonds(atommap *, atommap *);
    int mapChiral(atommap *, atommap *);
    int mapByIndex(atommap *, atommap *);
    int mapUniqueRefToTgt(atommap *, atommap *, int);

    int *AMap;
    bool maponly;
    Frame *newFrame;
    AmberParm *newParm;
    AmberParm *stripParm; // For stripping reference
  public:
    AtomMap(); 
    ~AtomMap();
    int init();
    int setup();
    int action();
};
#endif
