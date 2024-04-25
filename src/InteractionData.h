#ifndef INC_INTERACTIONDATA_H
#define INC_INTERACTIONDATA_H
#include "DataSet.h" // DataType
#include <map>
class DataFile;
class DataSetList;
class MetaData;
namespace Cpptraj {
/// Hold data sets corresponding to interactions between entities.
class InteractionData {
  public:
    InteractionData();

    /// Add interaction pair data set
    DataSet* AddInteractionSet(DataSetList&, DataSet::DataType, MetaData const&, int, int, DataFile*);
  private:
    typedef std::pair<int,int> Ipair;
    typedef std::pair<Ipair,DataSet*> PairType;
    typedef std::map<Ipair,DataSet*> MapType;

    MapType setMap_; ///< Hold map of interaction #s to data sets
};
}
#endif
