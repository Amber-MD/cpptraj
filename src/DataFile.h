#ifndef INC_DATAFILE_H
#define INC_DATAFILE_H
#include "DataIO.h"
#include "FileTypes.h"
/// Write DataSets to a file with specific format.
class DataFile {
    /// Allocator and description for file types. 
    static const FileTypes::AllocToken DF_AllocArray[];
    /// For associating keywords/extensions with file types. 
    static const FileTypes::KeyToken DF_KeyArray[];
    /// Keywords/extensions for types that support writes.
    static const FileTypes::KeyToken DF_WriteKeyArray[];
  public:
    /// Known data file formats.
    enum DataFormatType {
      DATAFILE=0, XMGRACE, GNUPLOT, XPLOR, OPENDX, REMLOG, MDOUT, EVECS,
      VECTRAJ, XVG, CCP4, CMATRIX, NCCMATRIX, CHARMMREPD, CHARMMFASTREP,
      CHARMMOUT, CPOUT, CHARMMRTFPRM, UNKNOWN_DATA 
    };
    DataFile();
    ~DataFile();
    // -------------------------------------------
    static void WriteHelp();
    /// List read options for each format.
    static void ReadOptions() { FileTypes::ReadOptions(DF_KeyArray,DF_AllocArray, UNKNOWN_DATA); }
    /// List write options for each format.
    static void WriteOptions(){ FileTypes::WriteOptions(DF_WriteKeyArray,DF_AllocArray,UNKNOWN_DATA); }
    /// \return Write format type from keyword in ArgList, or default
    static DataFormatType WriteFormatFromArg(ArgList& a, DataFormatType def) {
      return (DataFormatType)FileTypes::GetFormatFromArg(DF_WriteKeyArray,a,def);
    }
    /// \return string corresponding to format.
    static const char* FormatString(DataFormatType t) {
      return FileTypes::FormatDescription(DF_AllocArray, t);
    }
    /// \return string corresponding to file current format.
    const char* FormatString() const { return FileTypes::FormatDescription(DF_AllocArray,dfType_);}
    // -------------------------------------------
    /// Set debug level.
    void SetDebug(int);
    /// Set precision for all DataSets in DataFile
    void SetDataFilePrecision(int, int);
    /// Read data from DataFile to DataSets.
    int ReadDataIn(FileName const&, ArgList const&, DataSetList&);
    /// Read data from DataFile to DataSets; optionally append index to set name.
    int ReadDataIn(FileName const&, ArgList const&, DataSetList&, int, int);
    /// Read data from specific type of DataFile
    int ReadDataOfType(FileName const&, DataFormatType, DataSetList&);
    /// Set up DataFile for writing with optional args.
    int SetupDatafile(FileName const&, ArgList&, int);
    /// Set up DataFile for writing with specific format.
    int SetupDatafile(FileName const&, ArgList&, DataFormatType, int);
    /// Set up DataFile for writing to STDOUT (DataIO_Std) with optional arguments
    int SetupStdout(ArgList&, int);
    /// Set up DataFile for writing to STDOUT
    int SetupStdout(int d) { ArgList tmp; return SetupStdout(tmp, d); }
    /// Set up DataFile for writing, no args.
    int SetupDatafile(FileName const& f, int d) { ArgList a; return SetupDatafile(f, a, d); }
    /// Add a previously set-up DataSet to DataFile.
    int AddDataSet(DataSet*);
    /// Remove a set from the DataFile.
    int RemoveDataSet(DataSet*);
    /// Process DataFile-related arguments
    int ProcessArgs(ArgList&);
    int ProcessArgs(std::string const&); // TODO: Determine where this is used
    /// Write data in DataSets to disk.
    void WriteDataOut();
    /// \return string listing the names of all DataSets in DataFile.
    std::string DataSetNames() const;
    /// \return DataFile file name.
    FileName const& DataFilename() const { return filename_; }
    /// Used by DataFileList, indicates DataFile needs to be written. 
    void SetDFLwrite(bool fIn)           { dflWrite_ = fIn;  }
    /// Specify whether ensemble member number extension should be used.
    void SetEnsExt(bool b)               { ensExt_ = b;      }
    /// \return True if DataFile needs to be written.
    bool DFLwrite()                const { return dflWrite_; }
    /// \return DataFile format type.
    DataFormatType Type()          const { return dfType_;   }
  private:
    static DataIO* DetectFormat(FileName const&, DataFormatType&);
    int WriteSetsToFile(FileName const&, DataSetList&);
    int WriteWithEnsExtension();
    int WriteNoEnsExtension();

    int debug_;
    int dimension_;            ///< The dimension of all sets in the DataFile.
    DataFormatType dfType_;    ///< Format to read/write data in DataFile.
    bool dflWrite_;            ///< True: write file when DataFileList::WriteAllDF called.
    bool setDataSetPrecision_; ///< True: set default precision of incoming DataSets.
    bool sortSets_;            ///< True: Sort sets before write.
    bool ensExt_;              ///< If true append ensemble member number to file
    int default_width_;        ///< Default width of data sets added to this file.
    int default_precision_;    ///< Default precision of data sets added to this file.
    DataSetList SetList_;      ///< Array of pointers to associated DataSets.
    DataIO* dataio_;           ///< DataIO object for this DataFormatType.
    FileName filename_;        ///< DataFile file name.
    struct DimStruct {
      std::string label_;
      double min_;
      double step_;
    };
    /// Hold defaults for X, Y, and Z DataSet dimensions.
    std::vector<DimStruct> defaultDim_;
    /// True if min for X/Y/Z dim has been set.
    std::vector<bool> minIsSet_;
};
#endif
