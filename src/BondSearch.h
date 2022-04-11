#ifndef INC_BONDSEARCH_H
#define INC_BONDSEARCH_H
#include <map>
#include "Atom.h"
class Topology;
class Frame;
class Box;
class Vec3;
/** Used to search for bonds using coordinates in a Frame. */
class BondSearch {
  public:
    enum Type { SEARCH_REGULAR = 0, SEARCH_PAIRLIST, SEARCH_GRID, SEARCH_NONE };
    BondSearch();
    /// Search for bonds in given frame using specified search type, offset, and debug level.
    int FindBonds(Topology&, Type, Frame const&, double, int);
  private:
    /// Create a bounding box around atoms in frame, set vector with observed minimums
    static Box CreateBoundingBox(Frame const&, Vec3&);
    /// Create a bounding box around atoms in frame
    static Box CreateBoundingBox(Frame const&);
    /// Find/cache bond length squared for element pair (with offset)
    double GetBondLengthSquared(Atom::AtomicElementType, Atom::AtomicElementType, double);
    /// Search for bonds within residues
    void BondsWithinResidues(Topology&, Frame const&, double);
    /// \return Max # of bonds to given element
    static int MaxBonds(Atom::AtomicElementType);
    /// Search for bonds within residues, then bonds between residues using a grid.
    int GridSearch(Topology&, Frame const&, double, int);
    /// Search for bonds within residues, then bonds between adjacent residues.
    int ByResidueSearch(Topology&, Frame const&, double, int);
    /// Search for bonds using a pair list.
    int PairListSearch(Topology&, Frame const&, double, int);

    typedef std::pair<Atom::AtomicElementType,Atom::AtomicElementType> EltPair;
    typedef std::pair<EltPair, double> MapPair;
    typedef std::map<EltPair, double> BondMap;
    BondMap bond2Map_; ///< Used to cache squared bond lengths based on element types
};
#endif
