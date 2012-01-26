#ifndef INC_DATAFILELIST_H
#define INC_DATAFILELIST_H
#include "DataFile.h"
#include "DataSetList.h"
#include "ArgList.h"
#include <list>
// Class: DataFileList
/// Holds a list of DataFile classes. 
// NOTE: Currently implemented as a STL::list to make it easy to delete 
//       DataFiles, but this is currently not done. If that ability is
//       not needed a vector or an array would be more efficient.
class DataFileList {
    std::list<DataFile*> fileList;
    std::list<DataFile*>::iterator it;
    int debug;
  public:
    DataFileList();
    ~DataFileList();

    std::list<ArgList> DF_Args;

    void SetDebug(int);
    DataFile *GetDataFile(char *);
    DataFile *Add(char *, DataSet *);
    void Info();
    void Write();
    void ProcessDataFileArgs(DataSetList *);
};
#endif
