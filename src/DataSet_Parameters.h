#ifndef INC_DATASET_PARAMETERS_H
#define INC_DATASET_PARAMETERS_H
#include "DataSet.h"
#include "ParameterTypes.h"
#include "AtomTypeArray.h"
#include "ParameterHolders.h"
/// Hold parameters for atom types, bonds, angles, etc 
class DataSet_Parameters : public DataSet {
  public:
    DataSet_Parameters();
    static DataSet* Alloc() { return (DataSet*)new DataSet_Parameters(); }
    // ----- DataSet functions -------------------
    size_t Size()                                    const;
    void Info()                                      const;
    int Allocate(SizeArray const&)                         { return 1; }
    void Add(size_t, const void*)                          { return; }
    void WriteBuffer(CpptrajFile&, SizeArray const&) const { return; }
    int Append(DataSet*)                                   { return 1; }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
#   endif
    size_t MemUsageInBytes() const;
    // -------------------------------------------
    AtomTypeArray& AT()  { return atomTypes_; }
    ParmHolder<BondParmType>& BP() { return bondParm_; }
    ParmHolder<AngleParmType>& AP() { return angleParm_; }
    ParmHolder<BondParmType>& UB() { return ubParm_; }
    ParmHolder<DihedralParmType>& DP() { return dihParm_; }
    ParmHolder<DihedralParmType>& IP() { return impParm_; }

    ParmHolder<BondParmType> const& BP() const { return bondParm_; }
    ParmHolder<AngleParmType> const& AP() const { return angleParm_; }
    ParmHolder<BondParmType> const& UB() const { return ubParm_; }
    ParmHolder<DihedralParmType> const& DP() const { return dihParm_; }
    ParmHolder<DihedralParmType> const& IP() const { return impParm_; }

    void Debug() const;
  private:
    AtomTypeArray atomTypes_;
    ParmHolder<BondParmType> bondParm_;
    ParmHolder<AngleParmType> angleParm_;
    ParmHolder<BondParmType> ubParm_;
    ParmHolder<DihedralParmType> dihParm_;
    ParmHolder<DihedralParmType> impParm_;
};
#endif
