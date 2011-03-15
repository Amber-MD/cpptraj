#ifndef INC_ATOMMAP_H
#define INC_ATOMMAP_H
// AtomMap
#include "Action.h"

// Max number of bonds any given atom can have
#define MAXBONDS 6
// Hold 1 char for atom and each atom it is bonded to (+null)
// 1 + maxbonds + 1
#define ATOMIDLENGTH 8
// Hold atomID and each atomID it is bonded to (+null)
// atomidlength * (1 + maxbonds) + 1
#define UNIQUELENGTH 57
// Class: atommap
class atommap {
    struct mapatom {
      int bond[MAXBONDS];        // Holds indices of other bonded atoms
      int nbond;                 // Number of bonds
      bool complete;             // true: This atom an all bonded atoms have been mapped
      bool isChiral;             // true: Atom is a chiral center
      char atomID[ATOMIDLENGTH]; // ID created from this atom name, then bonded atom names 
      char unique[UNIQUELENGTH]; // ID created from this atomID, then bonded atomIDs
      bool isUnique;             // true: no other unique ID matches this atom
      bool isMapped;             // true: this atom has been mapped
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
    const char *atomID(int);
    const char *Aname(int);
    int setup();                    // Set up atom map, get atom elmt names, init map
    double getCut(char *, char *);  // Determine bond length cutoff between two atoms
    int calcDist();                 // Determine bonds between atoms based on distance
    void determineAtomID();         // Give each atom an id based on what atoms are bonded to it
    void markAtomComplete(int,bool);// Mark complete if atom and all atoms bonded to it are unique
    void markComplete();            // Mark unique atoms bonded to all unique atoms as complete
    bool BondIsRepeated(int,int);   // True if bonded atom has same name as another bonded atom
    // DEBUG
    void WriteMol2(char *);
};
// Class: AtomMap
class AtomMap : public Action {
    atommap RefMap;
    atommap TargetMap;
    int mapBondsToUnique(atommap *, atommap *);
    int mapChiral(atommap *, atommap *);
    int mapByIndex(atommap *, atommap *);
    int mapUniqueRefToTgt(atommap *, atommap *, int);
    int MapAtoms(atommap *, atommap *);

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
