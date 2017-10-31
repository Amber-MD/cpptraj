#ifndef INC_REPLICAINFO_H
#define INC_REPLICAINFO_H
#include <map>
/// Hold types/classes/functions having to do with processing ensemble/replica data.
namespace ReplicaInfo { // TODO just Replica
/// Replica sort target type
enum TargetType { NONE = 0, TEMP, INDICES, CRDIDX };

/// Hold temperature/indices for replica ensemble.
template <class T> class Map {
    typedef std::map<T, int> RmapType;
  public:
    Map() {}
    // Create a map from T array to replicas 0->N
    int CreateMap(std::vector<T> const&);
    T const& Duplicate() const { return duplicate_; }
    // Given T, find index in map
    int FindIndex( T const& ) const;
    typedef typename RmapType::const_iterator const_iterator;
    const_iterator begin() const { return repMap_.begin(); }
    const_iterator end()   const { return repMap_.end();   }
    bool empty()           const { return repMap_.empty(); }
    void ClearMap()              { repMap_.clear();        }
  private:
    RmapType repMap_;
    T duplicate_;
};
// Map::CreateMap()
template <class T> int Map<T>::CreateMap(std::vector<T> const& Vals) {
  std::set<T> tList;
  for (typename std::vector<T>::const_iterator val = Vals.begin(); val != Vals.end(); ++val)
  {
    std::pair<typename std::set<T>::iterator, bool> ret = tList.insert( *val );
    if (!ret.second) { // Duplicate value detected.
      duplicate_ = *val;
      return 1;
    }
  }
  repMap_.clear();
  // Values are now sorted lowest to highest in set.
  int repnum = 0;
  for (typename std::set<T>::const_iterator v0 = tList.begin(); v0 != tList.end(); ++v0, ++repnum)
    repMap_.insert(std::pair<T, int>(*v0, repnum));
  return 0;
}
// Map::FindIndex()
template <class T> int Map<T>::FindIndex( T const& Val ) const {
  typename RmapType::const_iterator rmap = repMap_.find( Val );
  if (rmap == repMap_.end())
    return -1;
  else
    return rmap->second;
}
} // END namespace Replica
#endif
