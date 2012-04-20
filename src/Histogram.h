#ifndef INC_HISTOGRAM_H
#define INC_HISTOGRAM_H
#include "Dimension.h"
#include "DataSetList.h"
// Class: Histogram
/// An N-dimensional histogram, the underlying data rep is 1D.
class Histogram {
  public:
    Histogram();
    /// Set debug level
    void SetDebug(int);
    /// Add dimension to histogram.  
    int AddDimension(Dimension&); 
    /// Bin given data. Dim must be the same as what the histogram has been set up for.
    int BinData(std::vector<double>&); 

    void PrintBins(char*, bool,bool);
    void Print_1D(DataSetList&);
    void Print_2D(DataSetList&);

    int CalcFreeE(double,int);
    int Normalize();

    /// Return given dimension
    Dimension& operator[](int);

    // Functions that return private vars
    int NumDimension() { return (int)dimensions_.size(); }
    int NumBins()      { return (int)Bins_.size();       }
  private:
    /// Debug level
    int debug_;
    /// Histogram dimensions
    std::vector<Dimension> dimensions_;
    /// Histogram data - double in case free E calculated
    std::vector<double> Bins_;
    /// 1 if data will wrap
    int isCircular_;
    /// Hold the bin indices in each dimension
    std::vector<int> BinIndices_;
    /// Set to true whenever highest order coord cycles in IncrementBinIndices
    bool hasCycled_;

    // Private Functions
    /// Set current bin and indices to 0
    void BinStart(bool);
    int BinIndicesToIndex();
    bool IncrementBinIndices();
};
#endif
