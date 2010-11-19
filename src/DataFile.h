#ifndef INC_DATAFILE_H
#define INC_DATAFILE_H
// DataFile
// Container for one or more DataSets that will be written out to a file.
#include "DataSet.h"
#include "PtrajFile.h"

class DataFile {
    int debug;
    DataSet **SetList; // Will point to addresses in DataSetList;
    int Nsets;
    bool noEmptyFrames;
    bool isInverted;  // If true sets will be written in rows instead of columns
    char *xlabel;

    void WriteGrace(PtrajFile *, int);
    void WriteData(PtrajFile *, int);
    void WriteDataInverted(PtrajFile *, int);
  public:
    char *filename;

    DataFile(char *);
    ~DataFile();

    void SetDebug(int);
    void SetXlabel(char*);
    void SetInverted();
    int AddSet(DataSet *);
    int NameIs(char *);
    void DataSetNames();
    void Write(int,bool);
};    

#endif
