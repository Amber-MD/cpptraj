#ifndef INC_DATASETLIST_H
#define INC_DATASETLIST_H
#include <vector>
#include "DataSet.h"
// Class: DataSetList
/// Hold list of data sets.
/** Main class for handling datasets. All dataset types can be allocated 
  * by DataSetList. DataSets are added to the list by various actions. 
  * There is a master DataSetList in Cpptraj that will hold data indexed 
  * by frame for use in analysis. All Data that may need to be available
  * for analysis should be in the master DataSetList.
  * This class is also used by DataFile to hold DataSets to be printed out.
  * This is done primarily to make use of the Get() functionality of 
  * DataSetList. In DataFile only copies of the sets are held, so they
  * are not freed by the destructor.
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
    /// Erase set from list
    void erase( const_iterator );
    /// Sort DataSets in list.
    void sort();
    /// True if no DataSets in list.
    bool empty()                  { return DataList_.empty();     }
    /// Return DataSet at didx.
    DataSet* operator[](int didx) { return DataList_[didx]; } // FIXME: No bounds check
    /// Return number of datasets in the list 
    int size()                    { return (int)DataList_.size(); }
    /// Return the max # expected frames
    int MaxFrames()               { return maxFrames_;            }
    /// Set DataSetList and underlying DataSet debug level
    void SetDebug(int);
    /// Set the max # frames expected to be read in. Used to preallocate DataSet size.
    void SetMax(int);
    /// Set width.precision of all DataSets in the list.
    void SetPrecisionOfDatasets(int, int);
    /// Separate input string into DataSet args.
    std::string ParseArgString(std::string const&, int&, std::string&);
    /// Get DataSet with specified name argument.
    DataSet* Get(const char *);
    /// Get DataSet with specified name, index, and aspect.
    DataSet* GetSet(std::string const&, int, std::string const&);
    /// Get multiple DataSets matching specified argument.
    DataSetList GetMultipleSets( std::string const& );
    /// Add DataSet to list with name, or default name if not specified. (TODO: OBSOLETE)
    DataSet* Add( DataSet::DataType, const char*, const char*);
    /// Add DataSet to list with name, or default name if not specified.
    DataSet* AddSet( DataSet::DataType, std::string const&, const char*);
    /// Add DataSet to list with name and index.
    DataSet* AddSetIdx( DataSet::DataType, std::string const&, int);
    /// Add DataSet to list with name and aspect.
    DataSet* AddSetAspect( DataSet::DataType, std::string const&, std::string const&);
    /// Add DataSet to list with name, idx, and aspect.
    DataSet* AddSetIdxAspect( DataSet::DataType, std::string const&, int, std::string const&);
    /// Add DataSet to list with name, index, aspect, and size.
    DataSet* AddSet( DataSet::DataType, std::string const&, int, std::string const&, int);
    /// Add DataSet to list that has already been allocated and Setup.
    int AddDataSet(DataSet*);
    /// Add a copy of the DataSet to the list; memory for DataSet will not be freed.
    void AddCopyOfSet(DataSet*);
    /// Print info on DataSets in the list
    void Info();
    /// Call sync for DataSets in the list (MPI only)
    void Sync();

    void VectorBegin();
    DataSet* NextVector();
    DataSet* NextMatrix();
    DataSet* NextModes();
  private:
    /// Dataset debug level
    int debug_;
    /// True if list contains copies that should not be freed in destructor.
    bool hasCopies_;

    typedef std::vector<DataSet*> DataListType;
    /// List of datasets
    DataListType DataList_;
    /// Expected number of frames to be read in.
    int maxFrames_;
    /// Used to iterate over vector datasets
    int vecidx_;

    struct dsl_cmp {
      inline bool operator()(DataSet* first, DataSet* second) const {
        return *first < *second;
      }
    };
};
#endif
