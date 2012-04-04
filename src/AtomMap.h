#ifndef INC_ATOMMAP_H
#define INC_ATOMMAP_H
#include "MapAtom.h"
#include "Topology.h"
class AtomMap {
  public:
    AtomMap();

    MapAtom& operator[](int);
    int Natom();
    void SetDebug(int);
    int Setup(Topology*);
    void ResetMapping();
    bool BondIsRepeated(int,int);
    void DetermineAtomIDs();
    void MarkAtomComplete(int);
    void CheckForCompleteAtoms();
    int CheckBonds();
    
  private:
    static MapAtom EMPTYMAPATOM;
    std::vector<MapAtom> mapatoms_;
    int debug_;
};
#endif
