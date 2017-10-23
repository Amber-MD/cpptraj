#ifndef INC_DATASET_PH
#define INC_DATASET_PH
#include "DataSet.h"
/// Hold data from constant pH simulations; protonation states of each residue.
class DataSet_PH : public DataSet {
    typedef std::vector<int> Iarray;
  public:
    DataSet_PH();
    static DataSet* Alloc() { return (DataSet*)new DataSet_PH(); }

    class Residue {
      public:
        Residue() {}
        Residue(int nframes) : states_(nframes, 0) {}
        void push_back(int state)     { states_.push_back(state); }
        typedef Iarray::const_iterator const_iterator;
        const_iterator begin() const { return states_.begin(); }
        const_iterator end()   const { return states_.end();   }
      private:
        Iarray states_;
    };
    typedef std::vector<Residue> Rarray;

    void AddState(unsigned int res, int state) {
      if (res >= residues_.size())
        residues_.resize(res+1, Residue(nframes_));
      residues_[res].push_back( state );
    }

    typedef Rarray::const_iterator const_iterator;
    const_iterator begin() const { return residues_.begin(); }
    const_iterator end()   const { return residues_.end();   }
    // -------------------------------------------
    size_t Size() const { return residues_.size(); }
    void Info()   const { return; }
    void WriteBuffer(CpptrajFile&, SizeArray const&) const { return; }
    int Allocate(SizeArray const&)  { return 0; } // TODO implement?
    void Add( size_t, const void* ) { return; }
    int Append(DataSet*)            { return 1; }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
#   endif
  private:

    Rarray residues_;
    unsigned int nframes_;
};
#endif
