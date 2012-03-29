#ifndef INC_DATASETLIST_H
#define INC_DATASETLIST_H
#include "DataSet.h"
#include <vector>
// Class: DataSetList
/// Hold list of data sets.
/** Main class for handling datasets. All dataset types can be allocated 
  * by DataSetList. Data sets are added to the list by various actions. 
  * There is a master DataSetList in Cpptraj that will hold data indexed 
  * by frame for use in analysis. Certain actions (e.g. NA structure analysis, 
  * closest etc) have datasets that are not indexed by frame; these actons will 
  * have their own separate dataset lists. By convention DataSetLists are the 
  * only structures that track DataSet memory; DataFiles should only hold the 
  * addresses of DataSets.
  */
class DataSetList {
  public:
    DataSetList();
    ~DataSetList();
    /// DataSetList default iterator
    typedef std::vector<DataSet*>::const_iterator const_iterator;
    /// Iterator to beginning of dataset list
    const_iterator begin() const;
    /// Iterator to end of dataset list
    const_iterator end() const;
    //size_t DataCharSize();
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
    DataSet *AddMulti(DataSet::DataType, char *, const char *);
    /// Add dataset to the list with name prefix_suffixN; set index.
    DataSet *AddMultiN(DataSet::DataType, const char *, const char *, int);
    char *checkName(char*, const char*);
    //DataSet *AddMatrix(char*, const char*, int, int);
    /// Add dataset to the list with given name
    DataSet *Add( DataSet::DataType, char*, const char*);
    DataSet *AddIdx( DataSet::DataType, char*, int);
    int AddDataSet(DataSet*);
    int PrepareForWrite();
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
  private:
    /// Dataset debug level
    int debug;
    /// If true, this data set list only points to data sets
    bool hasCopies_;
    /// List of datasets
    std::vector<DataSet*> DataList;
    /// Number of datasets
    int Ndata;
    /// Expected number of frames to be read in.
    int maxFrames;
    //int currentSet;
};
#endif
