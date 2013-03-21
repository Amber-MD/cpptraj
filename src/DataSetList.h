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
    void Clear();
    DataSetList& operator+=(DataSetList const&);
    /// DataSetList default iterator
    typedef std::vector<DataSet*>::const_iterator const_iterator;
    /// Iterator to beginning of dataset list
    const_iterator begin() const { return DataList_.begin();     }
    /// Iterator to end of dataset list
    const_iterator end()   const { return DataList_.end();       }
    /// True if no DataSets in list.
    bool empty()           const { return DataList_.empty();     }
    /// Return number of datasets in the list 
    int size()             const { return (int)DataList_.size(); }
    /// Return the max # expected frames
    int MaxFrames()        const { return maxFrames_;            }
    /// Erase set from list
    void erase( const_iterator );
    /// Sort DataSets in list.
    void sort();
    /// Return DataSet at didx.
    DataSet* operator[](int didx) { return DataList_[didx]; } // FIXME: No bounds check
    /// Set DataSetList and underlying DataSet debug level
    void SetDebug(int);
    /// Set maximum # elements expected to be read in.
    void SetMax(int);
    /// Allocate DataSet memory based on current maximum.
    void AllocateSets();
    /// Set width.precision of all DataSets in the list.
    void SetPrecisionOfDatasets(int, int);
    /// Get DataSet with specified name, index, and aspect.
    DataSet* GetSet(std::string const&, int, std::string const&) const;
    /// Get DataSet matching specified argument.
    DataSet* GetDataSet( std::string const& ) const;
    /// Get multiple DataSets matching specified argument.
    DataSetList GetMultipleSets( std::string const& ) const;
    /// Generate name based on given default and # of DataSets.
    std::string GenerateDefaultName(const char*) const;
    /// Add DataSet to list with name, or default name if not specified.
    DataSet* AddSet( DataSet::DataType, std::string const&, const char*);
    /// Add DataSet to list with name and index.
    DataSet* AddSetIdx( DataSet::DataType, std::string const&, int);
    /// Add DataSet to list with name and aspect.
    DataSet* AddSetAspect( DataSet::DataType, std::string const&, std::string const&);
    /// Add DataSet to list with name, idx, and aspect.
    DataSet* AddSetIdxAspect( DataSet::DataType, std::string const&, int, std::string const&);
    /// Add DataSet to list with name, idx, aspect, and legend.
    DataSet* AddSetIdxAspect( DataSet::DataType, std::string const&, int, std::string const&,
                              std::string const&);
    /// Add a copy of the DataSet to the list; memory for DataSet will not be freed.
    void AddCopyOfSet(DataSet*);
    /// Print info on DataSets in the list
    void List() const;
    /// Call sync for DataSets in the list (MPI only)
    void Sync();
    /// Find next set of specified type with given name.
    DataSet* FindSetOfType(std::string const&, DataSet::DataType) const;
    /// Find COORDS DataSet or create default COORDS DataSet.
    DataSet* FindCoordsSet(std::string const&);
  private:
    /// Separate input string into DataSet args.
    static std::string ParseArgString(std::string const&, std::string&, std::string&);

    typedef std::vector<DataSet*> DataListType;
    /// DataSet debug level
    int debug_;
    /// True if list contains copies that should not be freed in destructor.
    bool hasCopies_;
    /// List of DataSets
    DataListType DataList_;
    /// Expected number of frames to be read in.
    int maxFrames_;

    /// Used to sort DataSets
    struct dsl_cmp {
      inline bool operator()(DataSet* first, DataSet* second) const {
        return *first < *second;
      }
    };
};
#endif
