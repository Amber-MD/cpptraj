#ifndef INC_DATAIO_MDOUT_H
#define INC_DATAIO_MDOUT_H
#include "DataIO.h"
/// Read energies from Amber MDOUT files.
class DataIO_Mdout : public DataIO {
  public:
    DataIO_Mdout() {}
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_Mdout(); }
    static void ReadHelp();
    int processReadArgs(ArgList&) { return 0; }
    int ReadData(FileName const&,DataSetList&,std::string const&);
    int processWriteArgs(ArgList&) { return 0; }
    int WriteData(FileName const&, DataSetList const&)   { return 1; }
    bool ID_DataFormat(CpptrajFile&);
  private:
    typedef std::vector<std::string> Sarray;
    typedef std::vector<double> Darray;
    enum FieldType { Etot= 0, EPtot, GMAX, BOND,
                     ANGLE, DIHED, VDWAALS, EEL, EGB, EPB, ECAVITY, EDISPER,
                     VDW14, EEL14, RESTRAINT, EAMBER, Density,
                     RMS, EKtot, ESURF, EAMD_BOOST, VOLUME, TEMP,
                     PRESS, DVDL, N_FIELDTYPES };
    static FieldType getEindex(Sarray const&);
    int GetAmberEterms(const char*, Darray&, std::vector<bool>&);
    static const char* Enames[];
};
#endif
