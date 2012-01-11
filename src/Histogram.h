#ifndef INC_HISTOGRAM_H
#define INC_HISTOGRAM_H
/// Class: Histogram
/// An N-dimensional histogram, the underlying data rep is 1D.
class Histogram {
    int debug;
    /// dimensionType holds information about a histogram dimension
    struct dimensionType {
      char label[64];
      double min;
      double max;
      double step;
      int bins;
      int offset;
    };
    dimensionType *Dimension;
    int numDimension; ///< Number of dimensions

    double *Bins;    ///< Histogram data - double in case free E calculated
    int numBins;     ///< Total number of bins
    int currentBin;  ///< The current bin
    int isCircular;  ///< 1 if data will wrap
    int* BinIndices; ///< Hold the bin indices in each dimension

    // Private Functions
    void count2coord(int *);
  public:
    Histogram();
    ~Histogram();

    void SetDebug(int);

    int AddDimension(char*,double,double,double,int); // Add dimension to histogram with given
                                                      // min, max, step, and bins.
    int AddDimension(char*,double,double,double); // Add dimension to histogram with given
                                                  // min, max, and step. Bins will be calcd.
    int BinData(double*); // Bin the given data. Dimension must be the same as what
                          // the histogram has been set up for.
    int CalcFreeE(double,int);
    int Normalize();
    void PrintBins(bool,bool);

    void BinStart(bool);            // Set current bin and indices to 0
    void CurrentBinCoord(double *); // Set coordinates for current bin
    double CurrentBinData();        // Return value of current bin
    int NextBin();                  // Increment bin and indices
    int NBins(int);                 // Return number of bins for given dimension
    double Step(int);               // Return step for the given dimension
    double Min(int);                // Return min for the given dimension
    double Max(int);                // Return max for the given dimension
    char *Label(int);               // Return label for the given dimension
    double BinTotal();              // Return sum of all bin values

    // Functions that return private vars
    int NumDimension() { return numDimension; }
    int NumBins()      { return numBins;      }
    int CurrentBin()   { return currentBin;   }
};
#endif
