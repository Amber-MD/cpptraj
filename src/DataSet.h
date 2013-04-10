#ifndef INC_DATASET_H
#define INC_DATASET_H
#include <vector>
#include "CpptrajFile.h"
// Class: DataSet
/// Base class that all DataSet types will inherit.
/** All atomic classes inheriting the DataSet class must implement 8 routines:
  * Xmax, Size, FrameIsEmpty, Add, WriteBuffer, Sync, Dval, and CurrentDval 
  * (the last 2 not needed for String).
  * Other types wishing to use DataFile output should at least implement the 
  * Size, Xmax, and WriteBuffer routines.
  * DataSets are given certain attributes to make DataSet selection easier; 
  * these are name, index, and aspect. Name is typically associated with the
  * action that creates the dataset, e.g. RMSD or distance. Index is used
  * when and action outputs subsets of data, e.g. with RMSD it is possible to 
  * output per-residue RMSD, where the DataSet index corresponds to the residue
  * number. Aspect is used to further subdivide output data type; e.g. with 
  * nucleic acid analysis each base pair (divided by index) has shear,
  * stagger etc calculated.
  */
class DataSet {
  public:
    /// Type of data stored in DataSet
    enum DataType {
      UNKNOWN_DATA=0, DOUBLE, STRING, INT, FLOAT, VECTOR, MATRIX, MODES, 
      HIST, TRIMATRIX, MATRIX2D, COORDS
    };
    /// Source of data stored in DataSet, used by Analysis_Statistics
    enum scalarMode {
      UNKNOWN_MODE=0, M_DISTANCE, M_ANGLE, M_TORSION, M_PUCKER, M_RMS
    };
    /// Type of DataSet, used by Analysis_Statistics
    enum scalarType {
      UNDEFINED=0, ALPHA, BETA, GAMMA, DELTA, EPSILON, ZETA, PUCKER, CHI, 
      H1P,         C2P,   PHI,  PSI,   PCHI,  HBOND,   NOE
    };
    // 0           DIH    DIH   DIH    DIH    DIH      DIH   PUCK    DIH
    // DIH         DIH    DIH   DIH    DIH    DIST     DIST

    DataSet();          // Constructor
    DataSet(DataType,int,int,int);
    virtual ~DataSet() {} // Destructor - virtual since this class is inherited

    // ----------===== Inheritable functions =====----------
    virtual int Allocate(int)       { return 0; }
    /// Return the largest X/frame value added to set. 
    virtual int Xmax()              { return 0; }
    /// Return the number of data elements stored in the set.
    virtual int Size()              { return 0; }
    /// Used to check if a frame in DataSet has data.
    virtual int FrameIsEmpty(int)   { return 1; }
    /// Add data to the DataSet.
    /** A pointer to the data is passed in as void - it is up to the 
      * inheriting class to cast it. The X value for the data is passed 
      * in as well. It is expected that each successive X value will
      * be greater than the preceeding one (does not need to be
      * consecutive however). 
      */
    virtual void Add( int, void * ) { return;   }
    /// Write data at frame to file 
    virtual void WriteBuffer(CpptrajFile&,int)      { return; }
    /// Write 2D data to file
    virtual void Write2D(CpptrajFile&,int,int)      { return; }
    /// Return size of all dimensions
    // NOTE: Currently only used for 2D output
    virtual void GetDimensions( std::vector<int>& ) { return; }
    /// Consolodate this DataSet across all threads (MPI only)
    virtual int Sync()              { return 0; }
    /// Return data from data set as double precision
    virtual double Dval(int)        { return 0; }
    /// Return last value written to data set as double precision
    virtual double CurrentDval()    { return 0; } // TODO: Obsolete
    /// Print DataSet information
    virtual void Info()             { return;   }
    // -----------------------------------------------------
    // -----===== Public functions =====-----
    /// Set output precision
    void SetPrecision(int,int);
    /// Set up DataSet with given name and size
    int SetupSet(std::string const&,int,std::string const&);
    /// True if set is empty. 
    bool Empty();
    /// Used to set the data and header format strings 
    int SetDataSetFormat(bool);
    /// DataSet output label.
    std::string const& Legend();
    /// Set DataSet legend.
    void SetLegend( std::string const& lIn ) { legend_ = lIn; }
    /// Set scalar mode
    void SetScalar( scalarMode mIn ) { scalarmode_ = mIn; }
    /// Set scalar mode and type
    void SetScalar( scalarMode, scalarType );
    /// Check if name and/or index and aspect match this dataset.
    bool Matches(std::string const&, int, std::string const&);

    // -----===== Functions that return private vars =====-----
    /// Return DataSet base name.
    std::string const& Name()   const { return name_; }
    /// Return DataSet index.
    int Idx()                   const { return idx_; }
    /// Return DataSet aspect.
    std::string const& Aspect() const { return aspect_; }
    /// Return DataSet type.
    DataType Type()             const { return dType_; }
    /// Return DataSet type name.
    const char* TypeName()      const { return SetStrings[dType_]; }
    /// Return DataSet dimension.
    int Dim()                   const { return dim_; }
    /// Size in characters necessary to write data from this set.
    int Width()                 const { return width_; }
    /// Return scalar mode
    scalarMode ScalarMode()     const { return scalarmode_; }
    /// Return scalar type
    scalarType ScalarType()     const { return scalartype_; }
    /// Comparison for sorting, name/aspect/idx
    bool operator<(const DataSet& rhs) const {
      if ( name_ == rhs.name_ ) {
        if ( aspect_ == rhs.aspect_ ) {
          return ( idx_ < rhs.idx_ );
        } else {
          return ( aspect_ < rhs.aspect_ );
        }
      } else {
        return ( name_ < rhs.name_ );
      }
    }
  protected:
    std::string name_;         ///< Name of the DataSet
    int idx_;                  ///< DataSet index
    std::string aspect_;       ///< DataSet aspect.
    std::string legend_;       ///< DataSet legend.
    DataType dType_;           ///< The DataSet type
    int dim_;                  ///< DataSet dimension; used to determine how to write output.
    int width_;                ///< The output width of a data element
    int precision_;            ///< The output precision of a data element (if applicable)
    std::string format_;       ///< Output format of data
    const char *data_format_;  ///< Used to avoid constant calls to c_str
    scalarMode scalarmode_;    ///< Source of data in DataSet.
    scalarType scalartype_;    ///< Specific type of data in DataSet (if any).
  private:
    static const char* SetStrings[];
    bool GoodCalcType();
};
#endif 
