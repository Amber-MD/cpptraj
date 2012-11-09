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

    void SetDebug(int);
    DataFile* GetDataFile(std::string const&);
    DataFile* AddDataFile(std::string const&, ArgList&);
    DataFile* AddDataFile(std::string const&);
    DataFile* AddSetToFile(std::string const&,  DataSet*);
    void List();
    void Write();
    int ProcessDataFileArgs(ArgList&);
  private:
    typedef std::vector<DataFile*>::iterator df_iterator;
    std::vector<DataFile*> fileList_;
    int debug_;
};
#endif
