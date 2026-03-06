#ifndef INC_DATASET_LEAPOPTS_H
#define INC_DATASET_LEAPOPTS_H
#include "DataSet.h"
#include "Parm/GB_Params.h"
/// Hold LEaP options read from 'sourcing' a leaprc file in DataIO_LeapRC 
class DataSet_LeapOpts : public DataSet {
  public:
    DataSet_LeapOpts();
    static DataSet* Alloc() { return (DataSet*)new DataSet_LeapOpts(); }
    // ----- DataSet functions -------------------
    size_t Size()                                    const { return 0; }
    void Info()                                      const { return; }
    int Allocate(SizeArray const&)                         { return 1; }
    void Add(size_t, const void*)                          { return; }
    void WriteBuffer(CpptrajFile&, SizeArray const&) const { return; }
    int Append(DataSet*)                                   { return 1; }
    size_t MemUsageInBytes()                         const { return 0; }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
#   endif
    // -------------------------------------------
    int SetGbRadii(std::string const&);
    int SetSCEE(double);
    int SetSCNB(double);
    int SetDipoleDampFactor(double);
    int SetIpol(int);
    int SetFlexibleWater(bool);
    int SetDeleteExtraPointAngles(bool);

    Cpptraj::Parm::GB_RadiiType PbRadii() const { return pbradii_; }
    double SCEE() const { return scee_; }
    double SCNB() const { return scnb_; }
    double DipoleDampFactor() const { return dipoleDampFactor_; }
    int IPOL() const { return ipol_; }
    bool FlexibleWater() const { return flexibleWater_; }
    bool DeleteExtraPointAngles() const { return deleteExtraPointAngles_; }
  private:
    Cpptraj::Parm::GB_RadiiType pbradii_;
    double scee_;
    double scnb_;
    double dipoleDampFactor_;
    int ipol_;
    bool flexibleWater_;
    bool deleteExtraPointAngles_;
};
#endif
