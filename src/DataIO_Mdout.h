#ifndef INC_DATAIO_MDOUT_H
#define INC_DATAIO_MDOUT_H
#include "DataIO.h"
/// Read energies from Amber MDOUT files.
class DataIO_Mdout : public DataIO {
  public:
    DataIO_Mdout() {}
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_Mdout(); }
    static void ReadHelp();
    int ReadData(std::string const&,ArgList&,DataSetList&,std::string const&);
    int processWriteArgs(ArgList&) { return 0; }
    int WriteData(std::string const&, DataSetList const&) { return 1; }
    int WriteData2D(std::string const&, DataSet const&) { return 1; }
    int WriteData3D(std::string const&, DataSet const&) { return 1; }
    bool ID_DataFormat(CpptrajFile&);
  private:
    enum FieldType { NSTEP = 0, Etot, EPtot, GMAX, BOND,
                     ANGLE, DIHED, VDWAALS, EEL, EGB,
                     VDW14, EEL14, RESTRAINT, EAMBER, Density,
                     RMS, EKtot, ESURF, EAMD_BOOST, N_FIELDTYPES };
    static FieldType getEindex(std::vector<std::string> const&);
    static const char* Enames[];
};
#endif
