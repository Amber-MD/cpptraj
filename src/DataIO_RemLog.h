#ifndef INC_DATAIO_REMLOG_H
#define INC_DATAIO_REMLOG_H
#include "DataIO.h"
#include "BufferedLine.h"
/// Read replica exchange log data.
class DataIO_RemLog : public DataIO {
  public:
    DataIO_RemLog();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_RemLog(); }
    static void ReadHelp();
    int ReadData(std::string const&,ArgList&,DataSetList&,std::string const&);
    int processWriteArgs(ArgList&) { return 0; }
    int WriteData(std::string const&, DataSetList const&) { return 1; }
    int WriteData2D(std::string const&, DataSetList const&) { return 1; }
    int WriteData3D(std::string const&, DataSetList const&) { return 1; }
    bool ID_DataFormat(CpptrajFile&) { return false; }
  private:
    int MremdRead(std::vector<std::string> const&, DataSetList&, std::string const&);
 
    enum ExchgType { UNKNOWN = 0, TREMD, HREMD, MREMD };
    int ReadRemlogHeader(BufferedLine&, ExchgType&);
    int ReadRemdDimFile(std::string const&);
    int debug_;
    class GroupReplica;
    typedef std::vector<GroupReplica> GroupArray;
    typedef std::vector<GroupArray> GroupDimType;
    std::vector<GroupDimType> GroupDims_;
};

class DataIO_RemLog::GroupReplica {
  public:
    GroupReplica() : l_partner_(-1), me_(-1), r_partner_(-1) {}
    GroupReplica(const GroupReplica& rhs) :
      l_partner_(rhs.l_partner_), me_(rhs.me_), r_partner_(rhs.r_partner_) {}
    GroupReplica(int l, int m, int r) : l_partner_(l), me_(m), r_partner_(r) {}
    int L_partner() const { return l_partner_; }
    int Me()        const { return me_;        }
    int R_partner() const { return r_partner_; }
  private:
    int l_partner_, me_, r_partner_;
}; 
#endif
