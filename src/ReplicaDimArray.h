#ifndef INC_REPLICADIMARRAY_H
#define INC_REPLICADIMARRAY_H
#include <vector>
class ReplicaDimArray {
  public:
    ReplicaDimArray() {}
    ReplicaDimArray(ReplicaDimArray const& rhs) : remDims_(rhs.remDims_) {}
    ReplicaDimArray& operator=(ReplicaDimArray const& rhs) {
      if (this == &rhs) return *this;
      remDims_ = rhs.remDims_;
      return *this;
    }
    // NOTE: The indices for TEMPERATURE(1), HAMILTONIAN(3), and PH(4)
    //       currently match what is defined in Amber for REMD (remd.F90)
    //       EXCEPT RXSGLD, which uses the TEMPERATURE framework in Amber.
    //       Care should be taken to keep these in sync.
    enum RemDimType { UNKNOWN=0, TEMPERATURE, PARTIAL, HAMILTONIAN, PH, REDOX, RXSGLD };
    RemDimType DimType(int idx) const { return remDims_[idx]; }
    int operator[](int idx) const { return (int)remDims_[idx];         }
    int Ndims()             const { return (int)remDims_.size();       }
    bool empty()            const { return remDims_.empty();           }
    void AddRemdDimension(int d)         { remDims_.push_back((RemDimType)d); }
    void AddRemdDimension(RemDimType d)  { remDims_.push_back(d);             }
    void ChangeRemdDim(int d, RemDimType t) { remDims_[d] = t; }
    void clear()                         { remDims_.clear();                  }
    static const char* dimType(RemDimType type) {
      switch (type) {
        case UNKNOWN:     return "Unknown";     // 0
        case TEMPERATURE: return "Temperature"; // 1
        case PARTIAL:     return "Partial";     // 2 (UNUSED?)
        case HAMILTONIAN: return "Hamiltonian"; // 3
        case PH:          return "pH";          // 4
        case REDOX:       return "RedOx";       // 5
        case RXSGLD:      return "RXSGLD";      // 6 FIXME placeholder for future traj type
      }
      return 0; // Sanity check, should never reach.
    }
    const char* Description(int idx) const {
      if (idx >= 0 && idx < (int)remDims_.size()) return dimType( remDims_[idx] ); 
      return 0; // Sanity check, should never reach.
    }
    bool operator!=(const ReplicaDimArray& rhs) const {
      if (remDims_.size() != rhs.remDims_.size()) return true;
      std::vector<RemDimType>::const_iterator d1 = rhs.remDims_.begin();
      for (std::vector<RemDimType>::const_iterator d0 = remDims_.begin();
                                                   d0 != remDims_.end(); ++d0)
        if ( *d0 != *(d1++) ) return true;
      return false;
    }
#   ifdef MPI
    void assign( unsigned int n, RemDimType t ) { remDims_.assign(n, t); }
    int* Ptr() { return (int*)&remDims_[0]; }
#   endif
    size_t DataSize() const { return remDims_.size() * sizeof(RemDimType); }
  private:
    std::vector<RemDimType> remDims_;
};
#endif
