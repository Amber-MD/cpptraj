#ifndef INC_DATAIO_AMBERFF_H
#define INC_DATAIO_AMBERFF_H
#include "DataIO.h"
class DataSet_LeapOpts;
/// <Enter description of DataIO_AmberFF here>
class DataIO_AmberFF : public DataIO {
  public:
    DataIO_AmberFF();
    static void ReadHelp();
    static void WriteHelp();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_AmberFF(); }
    int processReadArgs(ArgList&);
    int ReadData(FileName const&, DataSetList&, std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&, DataSetList const&);
    bool ID_DataFormat(CpptrajFile&);
    /// Read in data with leap options
    int ReadData(FileName const&, DataSetList&, std::string const&, DataSet_LeapOpts*);
  private:

    std::string nbsetname_; ///< Nonbonded parameter set name to use
};
#endif
