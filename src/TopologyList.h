#ifndef INC_TOPOLOGYLIST_H
#define INC_TOPOLOGYLIST_H
#include "DispatchObject.h"
#include "Topology.h"
#include "FileList.h"
#include "ArgList.h"
// Class: TopologyList
/// Holds a list of Topology classes.
/** Can either add new topology by filename, or add existing topology by 
  * address. Can search for topology in list by index, full/base filename,
  * or tag.
  */
class TopologyList : public FileList {
  public:
    TopologyList();
    ~TopologyList();
    void Clear();
    Topology* GetParm(int) const;
    Topology* GetParm(ArgList&) const;
    int AddParmFile(std::string const&);
    int AddParmFile(std::string const&,std::string const&,bool,double);
    int AddParm(Topology*);
    void ReplaceParm(int, Topology*);
    void List() const;
  private:
    std::vector<Topology*> TopList_;
    bool hasCopies_;  ///< true: List contains addresses of topologies, do not delete
};
#endif
