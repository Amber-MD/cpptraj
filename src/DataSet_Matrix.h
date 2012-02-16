#ifndef INC_DATASET_MATRIX_H
#define INC_DATASET_MATRIX_H
#include <map>
#include "DataSet.h"
// Class: DataSet_Matrix
class DataSet_Matrix : public DataSet {
    std::map<int,double*> matrices;
    int n_rows;
    int m_cols;
    int msize;
    size_t MSIZE_BYTES;
  public:
    DataSet_Matrix();
    DataSet_Matrix(int,int);
    ~DataSet_Matrix();
    DataSet_Matrix(const DataSet_Matrix&);
    DataSet_Matrix &operator=(const DataSet_Matrix&);

    void Add(int, void*);
    int Get(void*, int);
};
#endif
