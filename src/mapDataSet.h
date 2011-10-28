#ifndef INC_MAPDATASET_H
#define INC_MAPDATASET_H
/// Class: mapDataSet
/// Hold a dataset which can have multiple columns, intended for 3 or 4 column
/// datasets, e.g. 2d histograms with or without stdev values.
#include <vector>
#include "DataSet.h"
class mapDataSet : public DataSet {
    std::vector<double> xData;
    std::vector<double> yData;
    std::vector<double> zData;
    int totalwidth; 
  public:
    mapDataSet();

    int Xmax();
    int isEmpty(int);
    void Add( int, void * );
    int Get(void *, int);
    char *Write(char *, int);
    void WriteBuffer(CharBuffer&,int);
    int Width();
    int Sync();
};
#endif
