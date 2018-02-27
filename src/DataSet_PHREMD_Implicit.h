#ifndef INC_DATASET_PH_IMPLICIT_H
#define INC_DATASET_PH_IMPLICIT_H
#include "DataSet_PHREMD.h"
/// Hold unsorted data from implicit constant pH REMD simulations 
class DataSet_PHREMD_Implicit : public DataSet_PHREMD {
  public:
    DataSet_PHREMD_Implicit();
    static DataSet* Alloc() { return (DataSet*)new DataSet_PHREMD_Implicit(); }
    // ----- DataSet functions -------------------
    size_t Size()                                    const { return records_.size(); }
    void Info()                                      const;
    int Allocate(SizeArray const&);
    void Add(size_t, const void*)                          { return; }
    void WriteBuffer(CpptrajFile&, SizeArray const&) const { return; }
    int Append(DataSet*)                                   { return 1; }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
#   endif
    // -------------------------------------------
    /// Hold single constant pH record
    class Record {
      public:
        Record() : pH_(-1.0), recType_(Cph::PARTIAL_RECORD) {}
        Record(float p, int r, Iarray const& s) : pH_(p), recType_(r) {
          if (recType_ >= 0)
            resStates_ = Iarray(1, s[r]);
          else
            resStates_ = s;
        }
      private:
        float pH_;         ///< solvent pH
        int recType_;      ///< Record type
        Iarray resStates_; ///< Hold state for all residues or given residue
    };
    typedef std::vector<Record> RecArray;
    /// \return Array of records.
    RecArray const& Records() const { return records_;       }
    /// Add record
    void AddRecord(Record const& s)  { records_.push_back(s); }
  private:
    RecArray records_;
};
#endif
