#ifndef INC_DATASET_PHREMD_H
#define INC_DATASET_PHREMD_H
#include "DataSet.h"
#include "Cph.h"
/// Base class for holding unsorted data from constant pH REMD simulations.
class DataSet_PHREMD : public DataSet {
  public:
    DataSet_PHREMD() {}
    DataSet_PHREMD(DataSet::DataType tIn, TextFormat const& fIn) :
      DataSet(tIn, PHREMD, fIn, 0) {} // 0 dim indicates DataSet-specific write
    virtual ~DataSet_PHREMD() {} // Virtual since inherited

    typedef std::vector<Cph::CpRes> Rarray;
    typedef Rarray::const_iterator const_iterator;
    /// Set residue information array
    void SetResidueInfo(Rarray const& r) { residues_ = r; }
    /// \return residue information array
    Rarray const& Residues()       const { return residues_; }
    /// \return specified residue info
    Cph::CpRes const& Res(int idx) const { return residues_[idx]; }
    /// Set Monte Carlo step size, initial time, and time step
    void SetTimeValues(Cph::CpTime const& t)       { time_ = t;    }
    /// \return constant pH time values
    Cph::CpTime const& Time()                const { return time_; }

    //void Resize(size_t); // TODO necessary?
  protected: // TODO private?
    typedef std::vector<int> Iarray;
    Rarray residues_;      ///< Array of residues.
    Cph::CpTime time_;     ///< Time values
};
#endif
