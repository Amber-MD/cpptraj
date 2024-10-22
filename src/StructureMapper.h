#ifndef INC_STRUCTUREMAPPER_H
#define INC_STRUCTUREMAPPER_H
#include "AtomMap.h"
class DataSet_Coords_REF;
/// Attempt to create a one-to-one mapping between atoms in two structures.
class StructureMapper {
  public:
    typedef std::vector<int> MapType;
    /// CONSTRUCTOR
    StructureMapper();
    /// DESTRUCTOR
    ~StructureMapper();

    /// Set the improper angle cutoff (in deg.) for mapping via chirality.
    void SetChiralImproperCut(double);

    int CreateMap(DataSet_Coords_REF*, DataSet_Coords_REF*, int);
    int CreateMapByResidue(DataSet_Coords_REF*, DataSet_Coords_REF*, int);
    int operator[](int idx)         const { return AMap_[idx];   }
    int MapSize()                   const { return (int)AMap_.size(); }
    MapType const& Map()            const { return AMap_;        }
    /// \return Number of target atoms successfully mapper to reference
    int Nmapped()                   const { return Nmapped_;     }
    /// \return true if all target atoms could be mapped.
    bool AllTgtMapped()             const { return (Nmapped_ == TgtMap_.Natom()); }
  private:
    int mapBondsToUnique(AtomMap&, AtomMap&);
    /// Map unmapped atoms bonded to mapped chiral cenmters, brute force
    int mapChiral_withUnmappedAtoms(MapType&, AtomMap&, AtomMap&,
                                    int, const int*, const int*,
                                    int, const int*, const int*) const;
    /// Get ref and tgt bonded atom priorities
    int getAtomPriorities(std::vector<int>&, std::vector<int>&, int, int);
    /// Map unmapped atoms bonded to mapped chiral centers via priority
    int mapChiral_viaPriority(MapType&, AtomMap&, AtomMap&, int, int);

    void map_atoms(int, int, AtomMap&, AtomMap&, int&, const char*);
    bool chiralImproperCutSatisfied(double) const;
    void check_large_delta(int, int, AtomMap const&, AtomMap const&, int, int, double) const;

    int mapChiral(AtomMap&, AtomMap&);
    int mapByIndex(AtomMap&, AtomMap&);
    int mapUniqueRefToTgt(AtomMap&, AtomMap&, int);
    int MapAtoms(AtomMap&, AtomMap&);
    int MapUniqueAtoms(AtomMap&, AtomMap&);
    int MapWithNoUniqueAtoms( AtomMap&, AtomMap& );
    void CountMappedAtoms();
    /// Clear pseudo topology/frames
    void clearPseudoTopFrame();

    AtomMap RefMap_; ///< Reference atom order.
    AtomMap TgtMap_; ///< Atoms to be reordered.
    MapType AMap_;
    int debug_;
    int Nmapped_; ///< Number of atoms in tgt actually mapped to ref
    Topology* refTop_; ///< Pseudo topology for current ref atoms being mapped to
    Topology* tgtTop_; ///< Pseudo topology for current tgt atoms being mapped
    Frame* refFrame_;  ///< Pseudo frame for current ref atoms being mapped to
    Frame* tgtFrame_;  ///< Pseudo frame for current tgt atoms being mapped
    double chiral_improper_cut_; ///< Cutoff for mapping chiral based on improper angle (radians)
};
#endif
