#ifndef INC_DATAIO_COORDS_H
#define INC_DATAIO_COORDS_H
#include "DataIO.h"
/// Load COORDS set from a file 
class DataIO_Coords : public DataIO {
  public:
    DataIO_Coords();
    static void ReadHelp();
    static void WriteHelp();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_Coords(); }
    int processReadArgs(ArgList&);
    int ReadData(FileName const&, DataSetList&, std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&, DataSetList const&);
    bool ID_DataFormat(CpptrajFile&);
  private:
    //bool is_parm_fmt_; ///< Set to true if format contains topology info
    //bool is_traj_fmt_; ///< Set to true if format contains coordinates
};
#endif
