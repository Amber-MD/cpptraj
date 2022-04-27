#ifndef INC_EXEC_COMPAREENERGY_H
#define INC_EXEC_COMPAREENERGY_H
#include "Exec.h"
#include "CharMask.h"
#include "OnlineVarT.h"
class DataSet_double;
/// <Enter description of Exec_CompareEnergy here>
class Exec_CompareEnergy : public Exec {
  public:
    Exec_CompareEnergy();
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_CompareEnergy(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    typedef std::vector<std::string> Sarray;

    static DataSet_Coords* GetCoordsSet(DataSetList const&, std::string const&);

    void BondEnergy(Frame const&, Topology const&,
                    Frame const&, Topology const&) const;
    int SetupBondArray(Topology const&, BondArray const&, BondArray const&);
    int SetupBondArrays(Topology const&, Topology const&);

    void AngleEnergy(Frame const&, Topology const&,
                    Frame const&, Topology const&) const;
    int SetupAngleArray(Topology const&, AngleArray const&, AngleArray const&);
    int SetupAngleArrays(Topology const&, Topology const&);

    int GetEnergies(DataSet_Coords*, DataSet_Coords*) const;

    CharMask mask1_;
    CharMask mask2_;
    CharMask mask3_;

    CpptrajFile* bondout_;
    DataSet_double* bondDeltaE_;
    DataSet_double* bondDeltaR_;
    BondArray commonBonds0_;
    BondArray commonBonds1_;
    Sarray bondNames_;

    CpptrajFile* angleout_;
    DataSet_double* angleDeltaE_;
    DataSet_double* angleDeltaR_;
    AngleArray commonAngles0_;
    AngleArray commonAngles1_;
    Sarray angleNames_;

};
#endif
