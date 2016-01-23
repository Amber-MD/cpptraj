#ifndef INC_DATASET_TOPOLOGY
#define INC_DATASET_TOPOLOGY
#include "DataSet.h"
#include "Topology.h"
/// Hold Topology data
class DataSet_Topology : public DataSet {
  public:
    DataSet_Topology() : DataSet(TOPOLOGY, GENERIC, TextFormat(), 0) {}
    static DataSet* Alloc() { return (DataSet*)new DataSet_Topology();}
    // ----- DataSet functions -------------------
    size_t Size()                  const { return (size_t)top_.Natom(); }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
#   endif
    void Info()                    const { top_.Brief(0);               }
    int Allocate(SizeArray const&)       { return 0;                    }
    void Add( size_t, const void* ) {}
    void WriteBuffer(CpptrajFile&, SizeArray const&) const {}
    int Append(DataSet*)                 { return 1;                    }
    // -------------------------------------------
    int LoadTopFromFile(ArgList const&, int);
    int StripTop( std::string const& );
    void SetTop(Topology const& t) { top_ = t;            }
    void SetPindex(int p)          { top_.SetPindex( p ); }
    Topology* TopPtr()             { return &top_; } // NOTE: pytraj currently relies on this 
    Topology const& Top() const    { return top_;  }
  private:
    Topology top_;
};
#endif
