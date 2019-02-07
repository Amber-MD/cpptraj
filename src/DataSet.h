#ifndef INC_DATASET_H
#define INC_DATASET_H
#include <cstddef> // size_t
#include <vector>
#include "MetaData.h"
#include "Dimension.h"
#include "AssociatedData.h"
#include "TextFormat.h"
#include "CpptrajFile.h"
#ifdef MPI
# include "Parallel.h"
#endif
/// Base class that all DataSet types will inherit.
/** The DataSet base class holds common information, like MetaData for
  * selection, TextFormat for text output, etc.
  * Note that specific parts of DataSet MetaData are not directly changeable
  * since changing things like name directly will affect searching. The
  * exception is ensemble number, which is not necessarily known at DataSet
  * creation time.
  */
class DataSet {
  public:
    typedef DataSet* (*AllocatorType)();
    typedef std::vector<size_t> SizeArray;
    /// DataSet base type. 
    enum DataType {
      UNKNOWN_DATA=0, DOUBLE, FLOAT, INTEGER, STRING, MATRIX_DBL, MATRIX_FLT, 
      COORDS, VECTOR, MODES, GRID_FLT, GRID_DBL, REMLOG, XYMESH, TRAJ, REF_FRAME,
      MAT3X3, TOPOLOGY, PH, PH_EXPL, PH_IMPL,
      PARAMETERS, PMATRIX_MEM, PMATRIX_NC
    };
    /// Group DataSet belongs to. TODO remove CLUSTERMATRIX
    enum DataGroup {
      GENERIC=0, SCALAR_1D, MATRIX_2D, GRID_3D, COORDINATES, CLUSTERMATRIX,
      PHREMD,    PWCACHE
    };

    DataSet();
    /// Set DataSet type, group, text output format, and # of dimensions.
    DataSet(DataType,DataGroup,TextFormat const&,int);
    DataSet(const DataSet&);
    DataSet& operator=(const DataSet&);
    /// Destructor. Virtual since this class is inherited.
    virtual ~DataSet();
    // ----------===== Inheritable functions =====----------
    /// \return the number of data elements stored in the DataSet.
    virtual size_t Size() const = 0;
    /// Print DataSet information //TODO return string instead?
    virtual void Info() const = 0;
    /// Write data to file given start indices. FIXME Buffer? Should this function take number of elements as well?
    virtual void WriteBuffer(CpptrajFile&, SizeArray const&) const = 0;
    /// \return value of coordinate for specified dimension d and position p.
    /** NOTE: It is assumed this can ALWAYS be represented as double precision. */
    virtual double Coord(unsigned int d, size_t p) const { return dim_[d].Coord(p); }
    /// Allocate data given numbers of elements.
    virtual int Allocate(SizeArray const&) = 0;
    /// Add element to data set.
    /** A pointer to the data is passed in as void - it is up to the
      * inheriting class to cast it. The X value for the data is passed
      * in as well. It is expected that each successive X value will
      * be greater than the preceeding one (does not need to be
      * consecutive however).
      */
    virtual void Add( size_t, const void* ) = 0;
    /// Can be used to append given data set to this one.
    virtual int Append(DataSet*) = 0;
    // TODO SizeInMB?
#   ifdef MPI
    /// Piece this DataSet together from multiple threads.
    virtual int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) = 0;
    // TODO pure virtual
    virtual int SendSet(int, Parallel::Comm const&) { return 1; }
    virtual int RecvSet(int, Parallel::Comm const&) { return 1; }
#   endif
    // -----------------------------------------------------
    /// Associate additional data with this set.
    void AssociateData(AssociatedData* a) { associatedData_.push_back( a->Copy() ); }
    /// Set DataSet MetaData
    int SetMeta(MetaData const&);
    /// Set DataSet ensemble number.
    void SetEnsemble(int e) { meta_.SetEnsembleNum( e ); }
    /// Set DataSet legend
    void SetLegend(std::string const& l) { meta_.SetLegend(l); }
    /// Set specific TextFormat part.
    TextFormat& SetupFormat() { return format_; }
    /// Set specified DataSet dimension.
    void SetDim(Dimension::DimIdxType i, Dimension const& d) { dim_[(int)i] = d; }
    void SetDim(int i, Dimension const& d)                   { dim_[i] = d;      }
#   ifdef MPI
    void SetNeedsSync(bool b) { needsSync_ = b;  }
    bool NeedsSync() const    { return needsSync_;  }
#   endif
    /// Check if name and/or index and aspect wildcard match this DataSet.
    bool Matches_WC(MetaData::SearchString const&, DataType) const;
    /// \return AssociateData of specified type.
    AssociatedData* GetAssociatedData(AssociatedData::AssociatedType) const;
    /// \return true if given metadata matches this set MetaData exactly.
    bool Matches_Exact(MetaData const& m) const { return meta_.Match_Exact(m); }
    /// \return True if DataSet is empty.
    bool Empty()                const { return (Size() == 0);      }
    /// \return DataSet output label.
    const char* legend()        const { return meta_.Legend().c_str();    }
    /// \return DataSet MetaData
    MetaData const& Meta()      const { return meta_; }
    /// \return DataSet TextFormat
    TextFormat const& Format()  const { return format_; }
    /// \return DataSet type.
    DataType Type()             const { return dType_;             }
    /// \return DataSet group
    DataGroup Group()           const { return dGroup_;            }

    /// \return number of dimensions.
    size_t Ndim()               const { return dim_.size();        }
    /// \return specified DataSet dimension. // TODO consolidate
    Dimension& ModifyDim(Dimension::DimIdxType i) { return dim_[(int)i]; }
    Dimension const& Dim(int i)             const { return dim_[i];      }

    /// Comparison for sorting, name/aspect/idx
    inline bool operator<(const DataSet& rhs) const { return meta_ < rhs.meta_; }
    /// Used to sort DataSet pointers (DataSetList, Array1D).
    struct DS_PtrCmp {
      inline bool operator()(DataSet const* first, DataSet const* second) const {
        return *first < *second;
      }
    };
  protected:
    TextFormat format_;         ///< Text output data format.
  private:
    /// Type to hold coordinate info for each dimension in DataSet.
    typedef std::vector<Dimension> DimArray;
    /// Type to hold any additional data associated with this data set.
    typedef std::vector<AssociatedData*> AdataArray;

    /// Clear any associated data.
    void ClearAssociatedData();
    // FIXME dim_ and associated functions like Coord need to be reworked
    //       depending on the set type. For example, dim_ doesnt really work
    //       for non-orthogonal grids.
    DimArray dim_;              ///< Holds info for each dimension in the DataSet.
    AdataArray associatedData_; ///< Holds any additonal data associated with this DataSet
    DataType dType_;            ///< The DataSet type
    DataGroup dGroup_;          ///< The DataSet group
    MetaData meta_;             ///< DataSet metadata
#   ifdef MPI
    bool needsSync_;            ///< True if DataSet needs sync. Should only be true once after run
#   endif
};
#endif 
