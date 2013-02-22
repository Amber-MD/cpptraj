#ifndef INC_ATOMMAP_H
#define INC_ATOMMAP_H
#include "MapAtom.h"
#include "Topology.h"
/// Used to set up mapping information for each atom.
class AtomMap {
  public:
    AtomMap();

    MapAtom& operator[](int);
    int Natom();
    void SetDebug(int);
    int Setup(Topology const&);
    int SetupResidue(Topology const&,int);
    void ResetMapping();
    bool BondIsRepeated(int,int);
    void DetermineAtomIDs();
    void MarkAtomComplete(int,bool);
    void CheckForCompleteAtoms();
    int CheckBonds();

  private:
    bool InvalidElement();

    static MapAtom EMPTYMAPATOM;
    std::vector<MapAtom> mapatoms_;
    int debug_;
};
#endif
