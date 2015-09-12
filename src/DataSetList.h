#ifndef INC_DATASETLIST_H
#define INC_DATASETLIST_H
#include <vector>
#include "DataSet.h"
#include "ArgList.h" // GetReferenceFrame
#include "ReferenceFrame.h" // GetReferenceFrame
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
    typedef std::vector<DataSet*> DataListType;
    typedef std::vector<double> Darray;

    DataSetList();
    ~DataSetList();

    void Clear();
    DataSetList& operator+=(DataSetList const&);
    /// \return DataSet at didx.
    DataSet* operator[](int didx) const { return DataList_[didx]; } // FIXME: No bounds check
    /// DataSetList default iterator
    typedef DataListType::const_iterator const_iterator;
    /// Iterator to beginning of dataset list
    const_iterator begin() const { return DataList_.begin(); }
    /// Iterator to end of dataset list
    const_iterator end()   const { return DataList_.end();   }
    /// True if no DataSets in list.
    bool empty()           const { return DataList_.empty(); }
    /// \return number of datasets in the list 
    size_t size()          const { return DataList_.size();  }
    /// \return Number of frames from last call to AllocateSets().
    long int MaxFrames()   const { return maxFrames_;        }
    /// Set current ensemble number.
    void SetEnsembleNum(int i)   { ensembleNum_ = i;         }
    /// Set DataSetList and underlying DataSet debug level
    void SetDebug(int);
    /// Set DataSets pending status.
    void SetDataSetsPending(bool b) { dataSetsPending_ = b; }
    /// Make all sets not part of an ensemble part of given ensemble.
    void MakeDataSetsEnsemble(int);
    /// \return Ensemble number; -1 if not an ensemble
    int EnsembleNum()      const { return ensembleNum_;      }
    /// \return True if Actions have indicated DataSets will be generated.
    bool DataSetsPending() const { return dataSetsPending_;  }

    /// Remove set from the list.
    void RemoveSet( DataSet* );
    /// Remove set from list but do not destroy.
    DataSet* PopSet( DataSet* );

    /// Allocate 1D DataSet memory based on current max# expected frames.
    void AllocateSets(long int);
    /// Set width and precision of specified DataSets in the list.
    void SetPrecisionOfDataSets(std::string const&, int, int);
 
   /// Get DataSet matching specified attributes.
    DataSet* CheckForSet( MetaData const& ) const;

    /// Get DataSet corresponding to specified argument.
    DataSet* GetDataSet( std::string const& ) const;
    /// Get multiple DataSets matching specified argument.
    DataSetList GetMultipleSets( std::string const& ) const;
    /// Select multiple sets, no warning if none found.
    DataSetList SelectSets( std::string const& ) const;
    /// Select multiple sets by group.
    DataSetList SelectGroupSets( std::string const&, DataSet::DataGroup ) const;

    /// Generate name based on given default and # of DataSets.
    std::string GenerateDefaultName(std::string const&) const;
    /// Add DataSet to list; set up default name if no name specified.
    DataSet* AddSet( DataSet::DataType, MetaData const&, const char*);
    /// Add DataSet to list with given MetaData.
    DataSet* AddSet( DataSet::DataType, MetaData const&);
    /// Add an already set up DataSet to list; memory for DataSet will be freed.
    int AddSet( DataSet* );
    /// Add new sets or append to existing ones.
    int AddOrAppendSets(Darray const&, DataListType const&);
    /// Add a copy of the DataSet to the list; memory for DataSet will not be freed.
    void AddCopyOfSet(DataSet*);

    /// Print info on DataSets in the list
    void List() const;
#   ifdef MPI
    /// Call sync for DataSets in the list (MPI only)
    void SynchronizeData();
#   endif

    /// Find next set of specified type with given name.
    DataSet* FindSetOfType(std::string const&, DataSet::DataType) const;
    /// Find COORDS DataSet or create default COORDS DataSet.
    DataSet* FindCoordsSet(std::string const&);
    /// reference arg help text
    static const char* RefArgs;
    /// Get reference frame DataSet from args
    ReferenceFrame GetReferenceFrame(ArgList&) const;
    /// List all reference frames.
    void ListReferenceFrames() const;
    /// topology arg help text
    static const char* TopArgs;
    // Get topology from args
    Topology* GetTopology(ArgList&) const;
    /// List all topologies
    void ListTopologies() const;
  private:
    DataSet* EraseSet( DataSet*, bool );
    /// Warn if DataSet not found but may be pending.
    inline void PendingWarning() const;
    /// Select sets according to argument and type.
    DataSetList SelectSets( std::string const&, DataSet::DataType ) const;
    /// Wrapper around DataList_.push_back() that does extra bookkeeping.
    void Push_Back(DataSet*);
    
    /// Hold number of frames from most recent AllocateSets() call.
    long int maxFrames_;
    /// DataSet debug level
    int debug_;
    /// Ensemble member number
    int ensembleNum_;
    /// True if list contains copies that should not be freed in destructor.
    bool hasCopies_;
    /// True if Actions will generate DataSets in the future.
    bool dataSetsPending_;
    /// List of DataSets
    DataListType DataList_;
    /// Pointers to reference data sets.
    DataListType RefList_;
    /// Pointers to topology data sets.
    DataListType TopList_;
    /// Hold descriptions and allocators for all DataSet types.
    struct DataToken {
      const char* Description;
      DataSet::AllocatorType Alloc;
    };
    static const DataToken DataArray[];
    typedef const DataToken* TokenPtr;
};
#endif
