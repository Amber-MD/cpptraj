#ifndef INC_DATAIO_MDOUT_H
#define INC_DATAIO_MDOUT_H
#include "DataIO.h"
#include <map>
/// Read energies from Amber MDOUT files.
class DataIO_Mdout : public DataIO {
  public:
    DataIO_Mdout();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_Mdout(); }
    static void ReadHelp();
    int processReadArgs(ArgList&) { return 0; }
    int ReadData(FileName const&,DataSetList&,std::string const&);
    int processWriteArgs(ArgList&) { return 0; }
    int WriteData(FileName const&, DataSetList const&)   { return 1; }
    bool ID_DataFormat(CpptrajFile&);
  private:
    typedef std::vector<double> Darray;
    typedef std::map<std::string, unsigned int> NameIdxMap;
    typedef std::pair<std::string, unsigned int> NameIdxPair;

    enum FieldType { ETOT= 0, EPTOT, GMAX, BOND,
                     ANGLE, DIHED, VDWAALS, EEL, EGB, EPB, ECAVITY, EDISPER,
                     VDW14, EEL14, RESTRAINT, EAMBER, DENSITY,
                     RMS, EKTOT, ESURF, EAMD_BOOST, VOLUME, TEMP,
                     PRESS, DVDL, N_FIELDTYPES };

    FieldType getTermIdx(std::string const&) const;
    int GetAmberEterms(const char*, Darray&, std::vector<bool>&);

    static const char* Enames_[];
    /// Map field names to indices into energy sets.
    NameIdxMap termIdxMap_;

};
#endif
