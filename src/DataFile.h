#ifndef INC_DATAFILE_H
#define INC_DATAFILE_H
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
    int ReadData(ArgList&, DataSetList&);
    int SetupDatafile(std::string const&, ArgList&, int);
    int AddSet(DataSet*);
    int ProcessArgs(ArgList&);
    int ProcessArgs(const char*);
    int ProcessArgs(std::string const&);
    void Write();

    void SetPrecision(int, int);
    FileName const& DataFilename() const  { return filename_; }
    void DataSetNames();
  private:
    int debug_;
    int dimension_;           ///< The dimension of all sets in the datafile.
    DataFormatType dataType_;
    bool isInverted_;
    DataSetList SetList_; 
    DataIO* dataio_;
    FileName filename_;

    void DetermineTypeFromExt( std::string const& );
    int SetupDataIO();
};
#endif
