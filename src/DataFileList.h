#ifndef INC_DATAFILELIST_H
#define INC_DATAFILELIST_H
#include "DispatchObject.h"
#include "FileList.h"
#include "DataFile.h"
#include "DataSet.h"
#include "ArgList.h"
// Class: DataFileList
/// Holds a list of DataFile classes. 
class DataFileList : public FileList {
  public:
    DataFileList();
    ~DataFileList();
    void Clear();
    void SetDebug(int);
    DataFile* GetDataFile(std::string const&) const;
    DataFile* AddDataFile(std::string const&, ArgList&);
    DataFile* AddDataFile(std::string const&);
    DataFile* AddSetToFile(std::string const&,  DataSet*);
    void List() const;
    void Write();
    int ProcessDataFileArgs(ArgList&);
  private:
    typedef std::vector<DataFile*> DFarray;
    DFarray fileList_;
    int debug_;
};
#endif
