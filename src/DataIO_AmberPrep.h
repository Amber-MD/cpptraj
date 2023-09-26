#ifndef INC_DATAIO_AMBERPREP_H
#define INC_DATAIO_AMBERPREP_H
#include "DataIO.h"
class BufferedLine;
/// Read in the Amber prep file format 
class DataIO_AmberPrep : public DataIO {
  public:
    DataIO_AmberPrep();
    static void ReadHelp();
    static void WriteHelp();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_AmberPrep(); }
    int processReadArgs(ArgList&);
    int ReadData(FileName const&, DataSetList&, std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&, DataSetList const&);
    bool ID_DataFormat(CpptrajFile&);
  private:
    int readAmberPrep(BufferedLine&, DataSetList&, std::string const&) const;

    bool removeDummyAtoms_; ///< Remove dummy atoms if true
};
#endif
