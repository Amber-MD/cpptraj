#ifndef INC_DATAFILELIST_H
#define INC_DATAFILELIST_H
#include "FileList.h"
#include "DataFile.h"
#include "DataSetList.h"
#include "ArgList.h"
// Class: DataFileList
/// Holds a list of DataFile classes. 
class DataFileList : public FileList {
  public:
    DataFileList();
    ~DataFileList();

    void SetDebug(int);
    DataFile* GetDataFile(std::string const&);
    DataFile* Add(const char*, DataSet*);
    DataFile* AddSetToFile(std::string const&,  DataSet*);
    void Info();
    void Write();
    void ProcessDataFileArgs(DataSetList*);
    void AddDatafileArg(ArgList&);
  private:
    typedef std::vector<DataFile*>::iterator df_iterator;
    std::vector<DataFile*> fileList_;
    std::vector<ArgList> DF_Args_;
    int debug_;
};
#endif
