#ifndef INC_DATASETLIST_H
#define INC_DATASETLIST_H
/// Class: DataSetList
/// Data sets are added to the list by various actions.
///   SetMax: Set the maximum number of frames expected to be read in. Used
///           to preallocate DataSet size.
///   Add: Add a dataset to the list with given type and name (use default 
///        name if given name is NULL).
///   Get: Return dataset with given name.
///   AddData: Add data to a dataset in the list with given index.
///   Info: Print names of datasets in the list.
///   Sync: Sync up datasets across all threads. 
#include "DataSet.h"
#include <vector>
class DataSetList {
    int debug;
    std::vector<DataSet*> DataList;
    int Ndata;
    int maxFrames;      // Expected number of frames that will be read in
    int currentSet;
  public:
    DataSetList();
    ~DataSetList();

    void SetDebug(int);
    void SetMax(int);
    DataSet *Get(char *);
    DataSet *Get(int);
    DataSet *GetDataSetN(int);
    DataSet *AddMulti(dataType, char *, const char *);
    DataSet *AddMultiN(dataType, const char *, const char *, int);
    DataSet *Add( dataType, char*, const char*);
    DataSet *AddIdx( dataType, char*, int);
    void Begin();
    int AddData(int, void *, int);
    int AddData(int, void *);
    //int AddDataToIdx(int, void *, int );
    void Info();
    void Sync();
    
    int Size() { return Ndata; }
};
#endif
