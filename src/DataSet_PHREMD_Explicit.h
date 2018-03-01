#ifndef INC_DATASET_PHREMD_EXPLICIT_H
#define INC_DATASET_PHREMD_EXPLICIT_H
#include "DataSet_PHREMD.h"
/// Hold unsorted data from explicit constant pH REMD simulations
class DataSet_PHREMD_Explicit : public DataSet_PHREMD {
    typedef std::vector<float> Farray;
  public:
    DataSet_PHREMD_Explicit();
    static DataSet* Alloc() { return (DataSet*)new DataSet_PHREMD_Explicit(); }
    // ----- DataSet functions -------------------
    size_t Size() const { return solvent_pH_.size(); }
    void Info()   const;
    void WriteBuffer(CpptrajFile&, SizeArray const&) const { return; }
    /// Reserve space for states of each residue
    int Allocate(SizeArray const&);
    void Add( size_t, const void* ) { return; }
    int Append(DataSet*)            { return 1; }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
    //void Reduce(Parallel::Comm const&, int);
#   endif
    // -------------------------------------------
    typedef Farray::const_iterator ph_iterator;
    /// \return Solvent pH array
    Farray const& pH_Values() const { return solvent_pH_; }
    /// \return Residue states array
    Iarray const& ResStates()         const { return resStates_;    }
    /// \return Record type for given index
    int RecordType(unsigned int idx)  const { return recType_[idx]; }
    /// Add residue states, pH and record type.
    void AddState(Iarray const& states, float pH, int recType) {
      for (Iarray::const_iterator it = states.begin(); it != states.end(); ++it)
        resStates_.push_back( *it );
      solvent_pH_.push_back( pH );
      recType_.push_back( recType );
    }
    /// \return First solvent pH value.
    float Initial_pH() const { return solvent_pH_.front(); }

    //void Resize(size_t); // TODO necessary?
  private:
    Farray solvent_pH_;    ///< Solvent pH values each frame
    Iarray recType_;       ///< Record type each frame.
    Iarray resStates_;     ///< State of each residue each frame: {R00, R10, R20}, {R01, ...}
};
#endif
