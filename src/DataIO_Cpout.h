#ifndef INC_DATAIO_CPOUT_H
#define INC_DATAIO_CPOUT_H
#include "DataIO.h"
#include "DataSet_PH.h"
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

    typedef std::vector<double> Darray;
    typedef std::vector<int> Iarray;

    struct StateInfo {
      int num_states_;
      int first_atom_;
      int num_atoms_;
      int first_state_;
      int first_charge_;
    };
    typedef std::vector<StateInfo> StateArray;

    int ReadCpin(FileName const&);

    FileName cpin_file_;
    FileType type_;
    int trescnt_;
    Darray charges_;
    DataSet_PH::Rarray Residues_;
};
#endif
