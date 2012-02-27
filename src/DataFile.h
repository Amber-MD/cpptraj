#ifndef INC_DATAFILE_H
#define INC_DATAFILE_H
#include <vector>
#include <string>
#include "DataSet.h"
#include "CpptrajFile.h"
// Class: DataFile
/// Container for one or more DataSets that will be written out to a file. 
/** Only the addresses of the data sets are stored in the datafile; the
  * actual datasets should reside somewhere else in memory. 
  */
class DataFile {
    int debug;
    std::vector<DataSet*> SetList;  ///< Will point to addresses in a DataSetList;
    int Nsets;
    bool noEmptyFrames;   ///< If true, frames in which no sets have data will be skipped
    bool isInverted;      ///< If true sets will be written in rows instead of columns
    bool noXcolumn;       ///< If true the Frame column will not be written.
    bool writeHeader;     ///< If true write header line for this data file.
    std::string xlabel;   ///< X axis label for grace plots
    int xcol_width;       ///< Width in chars of the X column
    int xcol_precision;   ///< Precision of the X column
    std::string x_format; ///< Format string for printing x coord
    std::string ylabel;   ///< Y axis label for grace plots
    int maxFrames;        ///< The largest X value of any sets in SetList
    std::string filename; ///< DataFile filename

    double xmin;
    double xstep;
    double ymin;
    double ystep;
    // The following are currently used for gnuplot output only
    bool useMap;
    bool printLabels;

    void SetupXcolumn();
    void WriteGrace(CpptrajFile *);
    void WriteGraceInverted(CpptrajFile *);
    void WriteData(CpptrajFile *);
    void WriteDataInverted(CpptrajFile *);
    void WriteGnuplot(CpptrajFile *);
  public:

    DataFile();
    DataFile(char *);

    void SetDebug(int);
    void SetXlabel(char*);
    void SetYlabel(char*);
    void SetInverted();
    void SetNoXcol();
    void SetNoHeader();
    void SetPrecision(char *, int, int);
    void SetNoEmptyFrames();
    void SetCoordMinStep(double,double,double,double);
    void SetXstep(double);
    void SetMap();
    void SetNoLabels();

    int AddSet(DataSet *);
    bool DataFileNameIs(char *);
    void DataSetNames();
    void Write();

    const char *Filename() { return filename.c_str(); }
};    
#endif
