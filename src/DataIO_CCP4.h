#ifndef INC_DATAIO_CCP4_H
#define INC_DATAIO_CCP4_H
#include "DataIO.h"
/// Write CCP4 format data files.
/** Developed based on dcoumentation found at:
  *   http://www.ccp4.ac.uk/html/maplib.html#description
  */
class DataIO_CCP4 : public DataIO {
  public:
    DataIO_CCP4() : DataIO(false, false, true) {} // Valid for 3D only
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_CCP4(); }

    bool ID_DataFormat(CpptrajFile&);
    int processReadArgs(ArgList&) { return 0; }
    int ReadData(FileName const&, DataSetList&, std::string const&);
    static void WriteHelp();
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&,DataSetList const&);
  private:
    union headerbyte { unsigned char c[224]; int i[56]; float f[56]; };

    static const size_t wSize; ///< Size of a word (currently 4 bytes)
    static bool MapCharsValid(const unsigned char*);
    int WriteSet3D(DataSetList::const_iterator const&, CpptrajFile&);
};
#endif
