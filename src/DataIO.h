#ifndef INC_DATAIO_H
#define INC_DATAIO_H
#include "ArgList.h"
#include "DataSetList.h"
#include "CpptrajFile.h"
// Class: DataIO
/// Base class that all DataIO objects inherit from.
class DataIO {
  public:
    virtual ~DataIO() {}
    /// Used to allocate DataIO object in DataFile.
    typedef DataIO* (*AllocatorType)();
    /// Used for help text.
    typedef void (*HelpType)();
    // ----- Inherited Functions -----------------
    virtual int ReadData(std::string const&,ArgList&,DataSetList&,std::string const&) = 0;
    virtual int processWriteArgs(ArgList&) = 0;
    virtual int WriteData(std::string const&, DataSetList const&) = 0;
    virtual int WriteData2D(std::string const&, DataSet const&) = 0;
    virtual int WriteData3D(std::string const&, DataSet const&) = 0;
    virtual int WriteDataInverted(std::string const&, DataSetList const &) = 0;
    virtual bool ID_DataFormat(CpptrajFile&) = 0;
  protected:
    // TODO: Move this to DataSet?
    static std::string SetupCoordFormat(size_t, Dimension const&, int, int);
    static Dimension DetermineXdim( std::vector<double> const& );
};
#endif
