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
  private:
    int init();
    int setup();
    int action();

    AtomMap RefMap_;
    Frame* RefFrame_;
    Topology* RefParm_;
    AtomMap TgtMap_;
    Frame* TgtFrame_;
    Topology* TgtParm_ ;

    std::vector<int> AMap_;
    bool maponly_;
    Frame* newFrame_;
    Topology* newParm_;
    Topology* stripParm_; // For stripping reference

    Frame rmsRefFrame_;
    Frame rmsTgtFrame_;
    bool rmsfit_;
    DataSet* rmsdata_;

    int mapBondsToUnique(AtomMap&, AtomMap&);
    int mapChiral(AtomMap&, AtomMap&);
    int mapByIndex(AtomMap&, AtomMap&);
    int mapUniqueRefToTgt(AtomMap&, AtomMap&, int);
    int MapAtoms(AtomMap&, AtomMap&);
    int MapUniqueAtoms(AtomMap&, AtomMap&);
    int MapWithNoUniqueAtoms( AtomMap&, AtomMap& );
};
#endif
