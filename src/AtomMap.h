#ifndef INC_ATOMMAP_H
#define INC_ATOMMAP_H
// AtomMap
#include "Action.h"

#define ATOMIDLENGTH 26
#define UNIQUELENGTH 131

typedef struct {
  int bond[4];      // Holds indices of other bonded atoms
  int nbond;        // Number of bonds
  int visited;      // 1=Has been visited by traceBonds
                    // ID created from this atom name, then bonded atom names alphabetically
  char atomID[ATOMIDLENGTH]; 
                    // ID created from this atomID, then bonded atomIDs
  char unique[UNIQUELENGTH];   
  int isUnique;     // 1 if no other unique ID matches this atom
} mapatom;

class atommap {
  public:
    mapatom *M;
    int natom;
    char **names;
    Frame *F;
    AmberParm *P;
    int debug;

    atommap();
    ~atommap();
    int setup();
    double getCut(char *, char *);
    int calcDist();
    void determineAtomID();
};

class AtomMap : public Action {
    atommap RefMap;
    atommap TargetMap;
  public:
    //AtomMap() {currentType=ATOMMAP;}
    AtomMap() {}
    ~AtomMap(){}
    int init();
};
#endif
