#ifndef INC_DATAFILELIST_H
#define INC_DATAFILELIST_H
// DataFileList
#include "DataFile.h"
#include <list>
//using namespace std;

class DataFileList : public std::list<DataFile*> {
    std::list<DataFile*>::iterator it;
    int debug;
  public:
    DataFileList();
    ~DataFileList();

    void SetDebug(int);
    DataFile *GetDataFile(char *);
    int Add(char *, DataSet *);
    void Info();
    void Write(int);
};
#endif
