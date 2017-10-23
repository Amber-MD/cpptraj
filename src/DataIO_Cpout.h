#ifndef INC_DATAIO_CPOUT_H
#define INC_DATAIO_CPOUT_H
#include "DataIO.h"
/// Read Amber cpout file 
class DataIO_Cpout : public DataIO {
  public:
    DataIO_Cpout();
    static void ReadHelp();
    static void WriteHelp();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_Cpout(); }
    int processReadArgs(ArgList&);
    int ReadData(FileName const&, DataSetList&, std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&, DataSetList const&);
    bool ID_DataFormat(CpptrajFile&);
  private:
    static const char* FMT_REDOX_;
    static const char* FMT_PH_;
    enum FileType { PH = 0, REDOX, NONE };

    int ReadCpin(FileName const&);

    FileName cpin_file_;
    FileType type_;
};
#endif
