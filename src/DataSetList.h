#ifndef INC_DATASETLIST_H
#define INC_DATASETLIST_H
#include "DataSet.h"
#include <vector>
// Class: DataSetList
/// Hold list of data sets.
/** Main class for handling datasets. All dataset types can be allocated 
  * by DataSetList. Data sets are added to the list by various actions. 
  * There is a master DataSetList in CpptrajState that will hold data indexed 
  * by frame for use in analysis. Certain actions (e.g. NA structure analysis, 
  * closest etc) have datasets that are not indexed by frame; these actons will 
  * have their own separate dataset lists. By convention DataSetLists are the 
  * only structures that track DataSet memory; DataFiles should only hold the 
  * addresses of DataSets.
  */
class DataSetList {
    /// Dataset debug level
    int debug;
    /// List of datasets
    std::vector<DataSet*> DataList;
    /// Number of datasets
    int Ndata;
    /// Expected number of frames to be read in.
    int maxFrames;
    //int currentSet;
  public:
    DataSetList();
    ~DataSetList();
    /// Set dataset debug level
    void SetDebug(int);
    /// Set the maximum # frames expected to be read in. Used to preallocate DataSet size.
    void SetMax(int);
    /// Set width.precision of all datasets in the list.
    void SetPrecisionOfDatasets(int, int);
    //DataSet & operator[](int);
    DataSet *Get(char *);
    /// Get dataset with given index (added with AddMultiN)
    DataSet *GetDataSetIdx(int);
    /// Get the dataset at the given position in the list.
    DataSet *GetDataSetN(int);
    /// Add dataset to the list with name prefix_suffix
    DataSet *AddMulti(dataType, char *, const char *);
    /// Add dataset to the list with name prefix_suffixN; set index.
    DataSet *AddMultiN(dataType, const char *, const char *, int);
    char *checkName(char*, const char*);
    DataSet *AddMatrix(char*, const char*, int, int);
    /// Add dataset to the list with given name
    DataSet *Add( dataType, char*, const char*);
    DataSet *AddIdx( dataType, char*, int);
    //void Begin();
    /// Add data to the given set in the list
    int AddData(int, void *, int);
    //int AddData(int, void *);
    //int AddDataToIdx(int, void *, int );
    /// Print info on datasets in the list
    void Info();
    /// Call sync for datasets in the list (MPI only)
    void Sync();
    /// Return number of datasets in the list 
    int Size() { return Ndata; }
    /// Return the max # expected frames
    int MaxFrames() { return maxFrames; }
};
#endif
