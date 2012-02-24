#ifndef INC_DATASET_H
#define INC_DATASET_H
#include <string>
#include "CharBuffer.h"
/// Type of data stored in DataSet
enum dataType {UNKNOWN_DATA, DOUBLE, STRING, INT, XYZ, FLOAT, MATRIX};
/// Set a printf-style format string
int SetFormatString(std::string&, dataType, int, int, bool);
// Class: DataSet
/// Base class that all dataset types will inherit.
/** All classes inheriting the DataSet class must implement 6 routines:
  * isEmpty, Add, Get, WriteBuffer, Width, and Sync.
  */
class DataSet {
  protected:
    std::string name;   ///< Name of the dataset
    int idx;            ///< Dataset index
    dataType dType;     ///< The dataset type
    int N;              ///< Number of data elements
    int current;        ///< The current data element
    int width;          ///< The output width of a data element
    int precision;      ///< The output precision of a data element (if applicable)
    int leadingSpace;   ///< 0 if leftAligned, 1 otherwise
    std::string format; ///< Output format of data
    const char *data_format;  ///< Used to avoid constant calls to c_str
    std::string header_format;///< Output format of DataSet name

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
    /// Return data from data set as double precision
    virtual double Dval(int)                  { return 0; }

    // Psuedo-iterator functions
    virtual void Begin()           { return;       }
    virtual bool NextValue()       { return false; }
    virtual double CurrentValue()  { return 0;     }

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
    /// Write the dataset name to character buffer
    void WriteNameToBuffer(CharBuffer &);
    /// Check if set has been written to, set format string.
    int CheckSet();
    /// Used to set the data and header format strings 
    int SetDataSetFormat(bool);
    /// Data set capacity
    int Capacity() { return N; }

    // -----===== Functions that return private vars =====-----
    /// Dataset name
    char *Name()           { return (char*)name.c_str();  }
    /// Set dataset index
    void SetIdx(int idxIn) { idx = idxIn;  }
    /// Return dataset index
    int Idx()              { return idx;   }
    /// Return dataset type
    dataType Type()        { return dType; }
};
#endif 
