#ifndef INC_DATASET_H
#define INC_DATASET_H
#include <string>
#include "CharBuffer.h"
// Class: DataSet
/// Base class that all dataset types will inherit.
/** All classes inheriting the DataSet class must implement 6 routines:
  * isEmpty, Add, Get, WriteBuffer, Width, and Sync.
  */
class DataSet {
  public:
    /// Base type of data stored in DataSet
    enum DataType {
      UNKNOWN_DATA=0, DOUBLE, STRING, INT, FLOAT, VECTOR, MATRIX, MODES
    };
    /// Source of data stored in DataSet
    enum scalarMode {
      UNKNOWN_MODE=0, M_DISTANCE, M_ANGLE, M_TORSION, M_PUCKER, M_RMS
    };
    /// Type of DataSet
    enum scalarType {
      UNDEFINED=0, ALPHA, BETA, GAMMA, DELTA, EPSILON, ZETA, PUCKER, CHI, 
      H1P,         C2P,   PHI,  PSI,   PCHI,  HBOND,   NOE
    };
    // 0           DIH    DIH   DIH    DIH    DIH      DIH   PUCK    DIH
    // DIH         DIH    DIH   DIH    DIH    DIST     DIST

    DataSet();          // Constructor
    virtual ~DataSet(); // Destructor - virtual since this class is inherited

    // -----===== Inheritable functions =====-----
    /// Return the largest X/frame value added to set. 
    /** By convention this should be the last value added.
      */
    virtual int Xmax()              { return 0; }
    /// Return the number of data elements stored in the set.
    virtual int Size()              { return 0; }
    /// Used to check if a frame in dataset has data.
    virtual int FrameIsEmpty(int)   { return 1; }
    /// Add data to the dataset.
    /** A pointer to the data is passed in  as void - it is up to the 
      * inheriting class to cast it. The X value for the data is passed 
      * in as well. 
      */
    virtual void Add( int, void * ) { return;   }
    /// Get data from the dataset
    virtual int Get( void *, int )  { return 1; }
    /// Write data at frame to character buffer
    virtual void WriteBuffer(CharBuffer&,int) { return;   }
    /// Size in characters necessary to write data from this set.
    virtual int Width()             { return 0; }
    /// Consolodate this dataset across all threads (MPI only)
    virtual int Sync()              { return 0; }
    /// Return data from data set as double precision
    virtual double Dval(int)       { return 0; }
    // Psuedo-iterator functions
    virtual void Begin()           { return;       }
    virtual bool NextValue()       { return false; }
    virtual double CurrentValue()  { return 0;     }
    // -------------------------------------------

    // Calculation routines for atomic types (DOUBLE, FLOAT, INT)
    /// Return the average/stdev of all values in the set
    double Avg(double *);
    /// Return the smallest value added to the set.
    double Min();
    /// Return the largest value added to the set
    double Max();

    // -----===== Public functions =====-----
    /// Set output precision
    void SetPrecision(int,int);
    /// Set up dataset with given name and size
    int Setup(const char*,int);
    /// Print dataset information
    void Info();
    /// Check if set has been written to.
    int CheckSet();
    /// Used to set the data and header format strings 
    int SetDataSetFormat(bool);
    /// Write the dataset name to character buffer
    void WriteNameToBuffer(CharBuffer &);

    // -----===== Functions that return private vars =====-----
    /// Dataset name
    std::string Name()     { return name_; }
    /// Printf-compatible name
    const char* c_str()    { return name_.c_str(); }
    /// Set dataset index
    void SetIdx(int idxIn) { idx_ = idxIn; }
    /// Return dataset index
    int Idx()              { return idx_; }
    /// Return dataset type
    DataType Type()        { return dType_; }
    /// Return scalar mode
    scalarMode ScalarMode() { return scalarmode_; }
    /// Return scalar type
    scalarType ScalarType() { return scalartype_; } 
  protected:
    std::string name_;         ///< Name of the dataset
    int idx_;                  ///< Dataset index
    DataType dType_;           ///< The dataset type
    int width_;                ///< The output width of a data element
    int precision_;            ///< The output precision of a data element (if applicable)
    int leadingSpace_;         ///< 0 if leftAligned, 1 otherwise
    std::string format_;       ///< Output format of data
    const char *data_format_;  ///< Used to avoid constant calls to c_str
    std::string header_format_;///< Output format of DataSet name
    scalarMode scalarmode_;    ///< Source of data in dataset.
    scalarType scalartype_;    ///< Specific type of data in dataset (if any).
};
#endif 
