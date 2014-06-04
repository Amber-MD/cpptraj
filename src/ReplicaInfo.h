#ifndef INC_REPLICAINFO_H
#define INC_REPLICAINFO_H
#include <map>
// FIXME: Consolidate with Trajin_Multi
namespace ReplicaInfo {
  enum TargetType { NONE = 0, TEMP, INDICES, CRDIDX };
}
/// Hold temperature/indices for replica ensemble.
template <class T> class ReplicaMap {
  public:
    ReplicaMap() {}
    int CreateMap(std::vector<T> const&);
  private:
    typedef std::map<T, int> RmapType;
    RmapType repMap_;
    T duplicate_;
};
// ReplicaMap::CreateMap()
template<class T> int ReplicaMap<T>::CreateMap(std::vector<T> const& Vals) {
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
#endif
