#ifndef INC_DATASET_XYZ_H
#define INC_DATASET_XYZ_H
/// Class: DataSet_XYZ
/// Hold a dataset which can have multiple columns, intended for 3 or 4 column
/// datasets, e.g. 2d histograms with or without stdev values.
#include <vector>
#include "DataSet.h"
class DataSet_XYZ : public DataSet {
    std::vector<double> xData;
    std::vector<double> yData;
    std::vector<double> zData;
    int totalwidth; 
  public:
    DataSet_XYZ();

    int Xmax();
    int isEmpty(int);
    void Add( int, void * );
    int Get(void *, int);
    void WriteBuffer(CharBuffer&,int);
    int Width();
    int Sync();
};
#endif
