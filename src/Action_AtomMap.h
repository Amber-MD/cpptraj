#ifndef INC_ATOMMAP_H
#define INC_ATOMMAP_H
#include "Action.h"
// Class: atommap
/// Used to identify sections of a molecule which can then be mapped to another molecule
class atommap {
    struct mapatom {
      std::vector<int> bond;     ///< Holds indices of other bonded atoms
      int nbond;                 ///< Total number of bonds (eventually obsolete)
      bool isChiral;             ///< true: Atom is a chiral center
      std::string atomID;        ///< ID created from this atom name, then bonded atom names
      std::string unique;        ///< ID created from this atomID, then bonded atomIDs
      bool isUnique;             ///< true: no other unique ID matches this atom
      int Nduplicated;           ///< If !isUnique, how many times is uniqueID duplicated
      bool isMapped;             ///< true: this atom has been mapped
      bool complete;             ///< true: This atom an all bonded atoms have been mapped
    };
    int debug;
    static std::string NO_ID;    ///< Static reference for when atom # out of bounds.
  public:
    mapatom *M;         ///< Array of atoms, contains bond info etc
    int natom;          ///< Number of atoms in the map
    char *names;        ///< 1 char Name of each atom
    Frame *mapFrame;    ///< Hold atom coords
    AmberParm *mapParm; ///< Hold corresponding parm

    atommap();
    ~atommap();
    void SetDebug(int);
    const std::string &atomID(int);
    const char *Aname(int);
    int setup();                    ///< Set up atom map, get atom elmt names, init map
    //double getCut(char *, char *);  // Determine bond length cutoff between two atoms
    int calcDist();                 ///< Determine bonds between atoms based on distance
    void determineAtomID();         ///< Give each atom an id based on what atoms are bonded to it
    void markAtomComplete(int,bool);///< Mark complete if atom and all atoms bonded to it are unique
    void markComplete();            ///< Mark unique atoms bonded to all unique atoms as complete
    bool BondIsRepeated(int,int);   ///< True if bonded atom has same name as another bonded atom
    void ResetMapping();            ///< Reset any existing map information
    // DEBUG
    void WriteMol2(char *);
};
// Class: AtomMap
/// Action used to map one moleucle to another using atommap
class AtomMap : public Action {
    atommap RefMap;
    atommap TargetMap;
    int mapBondsToUnique(atommap *, atommap *);
    int mapChiral(atommap *, atommap *);
    int mapByIndex(atommap *, atommap *);
    int mapUniqueRefToTgt(atommap *, atommap *, int);
    int MapAtoms(atommap *, atommap *);
    int MapUniqueAtoms(atommap *, atommap *);
    int MapWithNoUniqueAtoms( atommap *, atommap * );

    int *AMap;
    bool maponly;
    Frame *newFrame;
    AmberParm *newParm;
    AmberParm *stripParm; // For stripping reference

    Frame rmsRefFrame;
    Frame rmsTgtFrame;
    bool rmsfit;
    DataSet *rmsdata;
  public:
    AtomMap(); 
    ~AtomMap();
    int init();
    int setup();
    int action();
};
#endif
