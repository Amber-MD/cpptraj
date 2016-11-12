#ifndef INC_STRUCTUREMAPPER_H
#define INC_STRUCTUREMAPPER_H
#include "AtomMap.h"
#include "DataSet_Coords_REF.h"
/// Attempt to create a one-to-one mapping between atoms in two structures.
class StructureMapper {
  public:
    typedef std::vector<int> MapType;
    StructureMapper() {}
    int CreateMap(DataSet_Coords_REF*, DataSet_Coords_REF*, int);
  private:
    int mapBondsToUnique(AtomMap&, AtomMap&);
    int mapChiral(AtomMap&, AtomMap&);
    int mapByIndex(AtomMap&, AtomMap&);
    int mapUniqueRefToTgt(AtomMap&, AtomMap&, int);
    int MapAtoms(AtomMap&, AtomMap&);
    int MapUniqueAtoms(AtomMap&, AtomMap&);
    int MapWithNoUniqueAtoms( AtomMap&, AtomMap& );

    AtomMap RefMap_; ///< Reference atom order.
    AtomMap TgtMap_; ///< Atoms to be reordered.
    MapType AMap_;
    Frame rmsRefFrame_;
    Frame rmsTgtFrame_;
    DataSet_Coords_REF* RefFrame_;
    DataSet_Coords_REF* TgtFrame_;
    int debug_;
};
#endif
