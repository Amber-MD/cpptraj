#ifndef INC_DATAIO_H
#define INC_DATAIO_H
#include "ArgList.h"
#include "DataSetList.h"
// Class: DataIO
/// Base class that all DataIO objects inherit from.
class DataIO {
  public:
    DataIO();
    virtual ~DataIO() { }
    DataIO(const DataIO&);
    DataIO &operator=(const DataIO&);

    void SetDebug(int);
    void SetMaxFrames(int);
 
    int ProcessCommonArgs(ArgList &); 

    virtual int ReadData(std::string const&, DataSetList &)          { return 1;}
    virtual int processWriteArgs(ArgList &)                          { return 0;}
    virtual int WriteData(std::string const&, DataSetList &)         { return 1;}
    virtual int WriteData2D(std::string const&, DataSet&)            { return 1;}
    virtual int WriteDataInverted(std::string const&, DataSetList &) { return 1;}
  protected:
    int maxFrames_;
    int debug_;
    bool hasXcolumn_;
    int xcol_width_;
    int xcol_precision_;
    std::string x_format_;
    std::string x_label_;
    double xmin_;
    double xstep_;
    double xoffset_;
    bool printEmptyFrames_;

    void SetupXcolumn();
    void PrintX(CpptrajFile& file, int frame) {
      file.Printf(x_format_.c_str(), (xstep_ * ((double)frame+xoffset_)) + xmin_);
    }
};
#endif
