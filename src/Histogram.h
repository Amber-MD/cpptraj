#ifndef INC_HISTOGRAM_H
#define INC_HISTOGRAM_H
/// Class: Histogram
/// An N-dimensional histogram, the underlying data rep is 1D.
class Histogram {
    int debug;
    // dimensionType holds information about a histogram dimension
    struct dimensionType {
      char label[64];
      double min;
      double max;
      double step;
      int bins;
      int offset;
    };
    dimensionType *Dimension;
    int numDimension; // Number of dimensions

    int *Bins; // Histogram data
    int numBins; // Total number of bins

    // Private Functions
    void count2coord(int *);
  public:
    Histogram();
    ~Histogram();

    void SetDebug(int);

    int AddDimension(char*,double,double,double,int); // Add dimension to histogram with given
                                                      // min, max, step, and bins.
    int BinData(double*); // Bin the given data. Dimension must be the same as what
                          // the histogram has been set up for.
    void PrintBins(int);
};
#endif
