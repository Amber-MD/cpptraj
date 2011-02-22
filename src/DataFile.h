#ifndef INC_DATAFILE_H
#define INC_DATAFILE_H
// DataFile
// Container for one or more DataSets that will be written out to a file.
#include "DataSet.h"
#include "PtrajFile.h"

class DataFile {
    int debug;
    DataSet **SetList;  // Will point to addresses in DataSetList;
    int Nsets;
    bool noEmptyFrames; // If true, frames in which no sets have data will be skipped
    bool isInverted;    // If true sets will be written in rows instead of columns
    bool noXcolumn;     // If true the Frame column will not be written.
    char *xlabel;       // X axis label for grace plots
    int maxFrames;      // The largest X value of any sets in SetList

    void WriteGrace(PtrajFile *);
    void WriteGraceInverted(PtrajFile *);
    void WriteData(PtrajFile *);
    void WriteDataInverted(PtrajFile *);
    void WriteGnuplot(PtrajFile *);
  public:
    char *filename;

    DataFile(char *);
    ~DataFile();

    void SetDebug(int);
    void SetXlabel(char*);
    void SetInverted();
    void SetNoXcol();
    int AddSet(DataSet *);
    int NameIs(char *);
    void DataSetNames();
    void Write(bool);
};    

#endif
