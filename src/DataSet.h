#ifndef INC_DATASET_H
#define INC_DATASET_H
#include "CharBuffer.h"
/// Type of data stored in DataSet
enum dataType {UNKNOWN_DATA, DOUBLE, STRING, INT, XYZ, FLOAT}; 
// Class: DataSet
/// Base class that all dataset types will inherit.
/** All classes inheriting the DataSet class must implement 7 routines:
  * Allocate, isEmpty, Add, Get, WriteBuffer, Width, and Sync.
  */
class DataSet {
  protected:
    char *name;        ///< Name of the dataset
    int idx;           ///< Dataset index
    dataType dType;    ///< The dataset type
    int N;             ///< Number of data elements
    int current;       ///< The current data element
    int width;         ///< The output width of a data element
    int precision;     ///< The output precision of a data element (if applicable)
    char *format;      ///< Format of output
    bool isDynamic;    /**< True : N is not known, reallocate as N increases
                            False: N is known, allocate for N */
    /// Used to initialize the dataset. 
    /** Allocate space for the number of frames to be read (passed in via
      * Setup()). If the value passed in is <= 0 (i.e. isDynamic) the dataset is
      * allocated dynamically. Allocate is called by Setup()
      */
      // NOTE: Currently all datasets make use of the C STL Map container, so
      //       no explicit allocation is required.
    virtual int Allocate( )      { return 0; }
    /// Used to set the output format string based on dType
    void setFormatString();

  public:
    DataSet();          // Constructor
    virtual ~DataSet(); // Destructor - virtual since this class is inherited
    // -----===== Inheritable functions =====-----
    /// Return the largest X/frame value added to set. 
    /** By convention this should be the last value added.
      */
    virtual int Xmax()                        { return 0; }
    /// Used to check if a frame in dataset has data.
    virtual int isEmpty(int)                  { return 0; }
    /// Add data to the dataset.
    /** A pointer to the data is passed in  as void - it is up to the 
      * inheriting class to cast it. The X value for the data is passed 
      * in as well. 
      */
    virtual void Add( int, void * )           { return;   }
    /// Get data from the dataset
    virtual int Get( void *, int )            { return 1; }
    /// Write data at frame to character buffer
    virtual void WriteBuffer(CharBuffer&,int) { return;   }
    /// Size in characters necessary to write data from this set.
    virtual int Width()                       { return 0; }
    /// Consolodate this dataset across all threads (MPI only)
    virtual int Sync()                        { return 0; }
    /// Return the smallest value added to the set.
    virtual double Min()                      { return 0; }
    /// Return the largest value added to the set
    virtual double Max()                      { return 0; }
    /// Return the average of all values in the set
    virtual double Avg()                      { return 0; }
    // -----===== Public functions =====-----
    /// Set output precision
    void SetPrecision(int,int);
    /// Set up dataset with given name and size
    int Setup(char*,int);
    /// Print dataset information
    void Info();
    /// Write the dataset name to character buffer
    void WriteNameToBuffer(CharBuffer &, bool);
    char *Name(char *,bool);
    int CheckSet();
    // -----===== Functions that return private vars =====-----
    char *Name()           { return name;  }
    void SetIdx(int idxIn) { idx = idxIn;  }
    int Idx()              { return idx;   }
    dataType Type()        { return dType; }
};
#endif 
