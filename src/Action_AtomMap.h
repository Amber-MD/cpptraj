#ifndef INC_ACTION_ATOMMAP_H
#define INC_ACTION_ATOMMAP_H
#include "Action.h"
#include "AtomMap.h"
// Class: Action_AtomMap
/// Action used to map one molecule to another using AtomMaps
class Action_AtomMap : public Action {
  public:
    Action_AtomMap(); 
    ~Action_AtomMap();
    int init();
    int setup();
    int action();
  private:
    AtomMap RefMap;
    Frame *RefFrame;
    Topology *RefParm;
    AtomMap TgtMap;
    Frame *TgtFrame;
    Topology *TgtParm;
    int mapBondsToUnique(AtomMap&, AtomMap&);
    int mapChiral(AtomMap&, AtomMap&);
    int mapByIndex(AtomMap&, AtomMap&);
    int mapUniqueRefToTgt(AtomMap&, AtomMap&, int);
    int MapAtoms(AtomMap&, AtomMap&);
    int MapUniqueAtoms(AtomMap&, AtomMap&);
    int MapWithNoUniqueAtoms( AtomMap&, AtomMap& );

    std::vector<int> AMap;
    bool maponly;
    Frame *newFrame;
    Topology *newParm;
    Topology *stripParm; // For stripping reference

    Frame rmsRefFrame;
    Frame rmsTgtFrame;
    bool rmsfit;
    DataSet *rmsdata;
};
#endif
