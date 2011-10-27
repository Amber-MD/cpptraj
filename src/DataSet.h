#ifndef INC_DATASET_H
#define INC_DATASET_H
/// Class: DataSet
/// The DataSet class is a base class that all dataset types will inherit.
/// All classes inheriting the DataSet class may implement up to 5 routines:
///   Allocate:
///     Used to initialize the dataset. The number of frames to be read
///     is passed in. If the value passed in is <= 0 the dataset is 
///     allocated dynamically. Allocate is called by Setup().
///     NOTE: Currently datasets make use of the C STL Map container, so
///     no explicit allocation is required.
///   isEmpty:
///     Used to check if a specific frame in the dataset has had any data 
///     added to it. If not, output may be skipped for the dataset if 
///     the noEmptyFrames option is specified.
///   Add:
///     Used to add data to the dataset. A pointer to the data is passed in
///     as void - it is up to the inheriting class to cast it. The X value
///     for the data is passed in as well. 
///   Write:
///     Format data at the specified frame and place into the character buffer
///     for output.
///   Width:
///     Report the size in characters necessary to write data from this set.
///   Sync :
///     MPI only - sum up this dataset across all threads.
///   Xmax:
///     Return the largest X/frame value added to set. By convention this should
///     be the last value added.
enum dataType {UNKNOWN_DATA, DOUBLE, STRING, INT, MAP, FLOAT}; 
class DataSet {
  protected:
    char *name;        // Name of the dataset
    int idx;           // Dataset index
    dataType dType;    // The dataset type
    int N;             // Number of data elements
    int current;       // The current data element
    int width;         // The output width of a data element
    int precision;     // The output precision of a data element (if applicable)
    char *format;      // Format of output
    bool isDynamic;    // isDynamic = True : N is not known, reallocate as N increases
                       //           = False: N is known, allocate for N
    // If not isDynamic, Allocate will reserve space for N data elements 
    virtual int Allocate( )      { return 0; }
    void setFormatString();

  public:
    DataSet();          // Constructor
    virtual ~DataSet(); // Destructor - virtual since this class is inherited
    // Inheritable functions
    virtual int Xmax()               { return 0; }
    virtual int isEmpty(int)         { return 0; }
    virtual void Add( int, void * )  { return;   }
    virtual int Get( void *, int )   { return 1; }
    virtual char *Write(char*, int)  { return 0; }
    virtual int Width()              { return 0; }
    virtual int Sync()               { return 0; }
    virtual double Min()             { return 0; }
    virtual double Max()             { return 0; }
    // Public functions
    void SetPrecision(int,int);
    int Setup(char*,int);
    void Info();
    char *Name(char *,bool);
    int CheckSet();
    // Functions that return private vars
    char *Name() { return name; }
    void SetIdx(int idxIn) { idx = idxIn; }
    int Idx() { return idx; }
    dataType Type() {return dType;}
};
#endif 
