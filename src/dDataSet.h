#ifndef INC_DDATASET_H
#define INC_DDATASET_H

#include "DataSet.h"

class dDataSet: public DataSet {
    double *Data;
    int *Xvalues;

  public:

    dDataSet();
    ~dDataSet();

    int Allocate();
    void Add( int, void * );
    //void Write(FILE*,int);
    char *Write(char*,int);
    int Sync();
};
#endif
