#ifndef INC_DATASET_H
#define INC_DATASET_H
/* DataSet
 * The DataSet class is a base class that all dataset types will inherit.
 * All classes inheriting the DataSet class may implement up to 5 routines:
 *   Allocate:
 *     Used to initialize the dataset. The number of frames to be read
 *     is passed in. If the value passed in is <= 0 the dataset is 
 *     allocated dynamically. Allocate is called by Setup().
 *     NOTE: Currently datasets make use of the C STL Map container, so
 *     no explicit allocation is required.
 *   isEmpty:
 *     Used to check if a specific frame in the dataset has had any data 
 *     added to it. If not, output may be skipped for the dataset if 
 *     the noEmptyFrames option is specified.
 *   Add:
 *     Used to add data to the dataset. A pointer to the data is passed in
 *     as void - it is up to the inheriting class to cast it. The X value
 *     for the data is passed in as well. 
 *   Write:
 *     Format data at the specified frame and place into the character buffer
 *     for output.
 *   Sync :
 *     MPI only - sum up this dataset across all threads. 
 */

enum dataType {
  UNKNOWN_DATA, DOUBLE, STRING, INT
};

class DataSet {
  protected:
    char *name;        // Name of the dataset
    int N;             // Number of data elements
    int current;       // The current data element
    int isDynamic;     // 1=N is not known, reallocate as N increases
                       // 0=N is known, allocate for N
    // Inherited by classes
    virtual int Allocate( )      { return 0; }

  public:

    DataSet();          // Constructor
    virtual ~DataSet(); // Destructor - virtual since this class is inherited

    virtual int isEmpty(int)         { return 0; }
    virtual void Add( int, void * )  { return;   }
    virtual char *Write(char*, int)   { return 0; }
    virtual int Sync()               { return 0; }

    int Setup(char*,int);
    void Info();
    char *Name() { return name; }
    int CheckSet();
};
#endif 
