#ifndef INC_MOL_H
#define INC_MOL_H
#include "Topology.h"
/// Routines associated with molecules.
namespace Mol {
  typedef std::vector<int> Iarray;
  /// Used to return unique mol types from UniqueCount
  class Type {
    public:
      Type(int i,int a,int r,std::string const& n) :
        idxs_(1,i), natom_(a), nres_(r), name_(n) {}
      void UpdateCount(int i) { idxs_.push_back(i); }
      Iarray idxs_;      ///< Molecule indices for molecules of this type
      int natom_;        ///< Number of atoms in molecule
      int nres_;         ///< Number of residues in molecule
      std::string name_; ///< Molecule type name
  };
  typedef std::vector<Type> Marray;
  Marray UniqueCount(Topology const&, CharMask const&);
  Marray UniqueCount(Topology const&);
}
#endif
