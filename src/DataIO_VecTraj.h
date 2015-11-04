#ifndef INC_DATAIO_VECTRAJ_H
#define INC_DATAIO_VECTRAJ_H
#include "DataIO.h"
#include "TrajectoryFile.h"
/// Write vector data to pseudo trajectory.
class DataIO_VecTraj : public DataIO {
  public:
    DataIO_VecTraj();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_VecTraj(); }
    static void WriteHelp();
    int processReadArgs(ArgList&) { return 0; }
    int ReadData(FileName const&,DataSetList&,std::string const&) {return 1;}
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&,DataSetList const&);
    bool ID_DataFormat(CpptrajFile&) { return false; }
  private:
    TrajectoryFile::TrajFormatType trajoutFmt_;
    std::string parmoutName_;
    bool includeOrigin_;
};
#endif
