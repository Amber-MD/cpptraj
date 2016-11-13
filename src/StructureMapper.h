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
    int operator[](int idx)         const { return AMap_[idx];   }
    int MapSize()                   const { return (int)AMap_.size(); }
    MapType const& Map()            const { return AMap_;        }
    MapAtom const& RefAtom(int idx) const { return RefMap_[idx]; }
    MapAtom const& TgtAtom(int idx) const { return TgtMap_[idx]; }
    /// \return Number of target atoms successfully mapper to reference
    int Nmapped()                   const { return Nmapped_;     }
    /// \return true if all target atoms could be mapped.
    bool AllTgtMapped()             const { return (Nmapped_ == TgtMap_.Natom()); }
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
    int debug_;
    int Nmapped_; ///< Number of atoms in tgt actually mapped to ref
};
#endif
