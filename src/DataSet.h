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
    /// Type of data stored in DataSet
    enum DataType {
      UNKNOWN_DATA, DOUBLE, STRING, INT, FLOAT
    };

    DataSet();          // Constructor
    virtual ~DataSet(); // Destructor - virtual since this class is inherited

    // -----===== Inheritable functions =====-----
    /// Return the largest X/frame value added to set. 
    /** By convention this should be the last value added.
      */
    virtual int Xmax()              { return 0; }
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
    int Setup(char*,int);
    /// Print dataset information
    void Info();
    /// Check if set has been written to.
    int CheckSet();
    /// Used to set the data and header format strings 
    int SetDataSetFormat(bool);
    /// Write the dataset name to character buffer
    void WriteNameToBuffer(CharBuffer &);

    // -----===== Functions that return private vars =====-----
    /// Data set capacity
    //int Capacity();
    /// Dataset name
    char *Name();
    /// Set dataset index
    void SetIdx(int);
    /// Return dataset index
    int Idx();
    /// Return dataset type
    DataType Type();
  protected:
    std::string name_;   ///< Name of the dataset
    int idx_;            ///< Dataset index
    DataType dType_;     ///< The dataset type
    //int N_;              ///< Number of data elements
    int current_;        ///< The current data element
    int width_;          ///< The output width of a data element
    int precision_;      ///< The output precision of a data element (if applicable)
    int leadingSpace_;   ///< 0 if leftAligned, 1 otherwise
    std::string format_; ///< Output format of data
    const char *data_format_;  ///< Used to avoid constant calls to c_str
    std::string header_format_;///< Output format of DataSet name
};
#endif 
