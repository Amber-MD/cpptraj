#ifndef INC_EXEC_COMPAREENERGY_H
#define INC_EXEC_COMPAREENERGY_H
#include "Exec.h"
#include "CharMask.h"
/// <Enter description of Exec_CompareEnergy here>
class Exec_CompareEnergy : public Exec {
  public:
    Exec_CompareEnergy();
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_CompareEnergy(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    static DataSet_Coords* GetCoordsSet(DataSetList const&, std::string const&);
    int GetEnergies(DataSet_Coords*, DataSet_Coords*) const;
    void CalcBondEnergy(Frame const&, BondArray const&, BondParmArray const&,
                        Frame const&, BondArray const&, BondParmArray const&) const;
    void BondEnergy(Frame const&, Topology const&,
                    Frame const&, Topology const&) const;

    CharMask mask1_;
    CharMask mask2_;
    CpptrajFile* bondout_;
};
#endif
