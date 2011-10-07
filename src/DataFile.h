#ifndef INC_DATAFILE_H
#define INC_DATAFILE_H
/// Class: DataFile
/// Container for one or more DataSets that will be written out to a file.
#include "DataSet.h"
#include "CpptrajFile.h"
class DataFile {
    int debug;
    DataSet **SetList;  // Will point to addresses in DataSetList;
    int Nsets;
    bool noEmptyFrames; // If true, frames in which no sets have data will be skipped
    bool isInverted;    // If true sets will be written in rows instead of columns
    bool noXcolumn;     // If true the Frame column will not be written.
    char *xlabel;       // X axis label for grace plots
    char *ylabel;       // Y axis label for grace plots
    int maxFrames;      // The largest X value of any sets in SetList

    // The following are currently used for gnuplot output only
    double xmin;
    double xstep;
    double ymin;
    double ystep;
    bool useMap;
    bool printLabels;

    void WriteGrace(CpptrajFile *);
    void WriteGraceInverted(CpptrajFile *);
    void WriteData(CpptrajFile *);
    void WriteDataInverted(CpptrajFile *);
    void WriteGnuplot(CpptrajFile *);
  public:
    char *filename;

    DataFile();
    DataFile(char *);
    ~DataFile();

    void SetDebug(int);
    void SetXlabel(char*);
    void SetYlabel(char*);
    void SetInverted();
    void SetNoXcol();
    void SetPrecision(char *, int, int);
    void SetNoEmptyFrames();
    void SetCoordMinStep(double,double,double,double);
    void SetMap();
    void SetNoLabels();

    int AddSet(DataSet *);
    int NameIs(char *);
    void DataSetNames();
    void Write();
};    
#endif
