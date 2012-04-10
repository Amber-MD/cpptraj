#ifndef INC_DATAFILE_H
#define INC_DATAFILE_H
#include "CpptrajFile.h"
#include "DataIO.h"
// Class: DataFile
// New datafile class for use with DataIO
/// Write datasets in a list to a file.
class DataFile {
  public:
    enum DataFormatType {
      UNKNOWN_DATA=0, DATAFILE, XMGRACE, GNUPLOT
    };

    DataFile();
    ~DataFile();

    void SetDebug(int);
    int SetupDatafile(char*);
    int AddSet(DataSet*);
    int ProcessArgs(ArgList&);
    int ProcessArgs(const char*);
    void SetXlabel(char*);
    void SetYlabel(char*);
    void SetCoordMinStep(double,double,double,double);
    void Write();

    void SetPrecision(char *, int, int);
    const char *Filename();
    void DataSetNames();
  private:
    int debug_;
    DataFormatType dataType_;
    bool isInverted_;
    DataSetList SetList_;
    DataIO *dataio_;
};
#endif
