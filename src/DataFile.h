#ifndef INC_DATAFILE_H
#define INC_DATAFILE_H
#include "DataIO.h"
#include "FileName.h"
// Class: DataFile
/// Write DataSets to a file with specific format.
class DataFile {
  public:
    /// Known data file formats.
    enum DataFormatType {
      DATAFILE, XMGRACE, GNUPLOT, XPLOR, OPENDX, REMLOG, UNKNOWN_DATA 
    };

    DataFile();
    ~DataFile();
    /// \return format type from keyword in ArgList.
    static DataFormatType GetFormatFromArg(ArgList&);
    /// \return format type from keyword.
    static DataFormatType GetFormatFromString(std::string const&);
    /// \return standard file extension for data file format.
    static std::string GetExtensionForType(DataFormatType);
    /// \return type from extension.
    static DataFormatType GetTypeFromExtension(std::string const&);

    void SetDebug(int);
    /// \return string corresponding to format.
    const char* FormatString() const;

    int ReadData(ArgList&, DataSetList&);
    int SetupDatafile(std::string const&, ArgList&, int);
    int AddSet(DataSet*);
    int ProcessArgs(ArgList&);
    int ProcessArgs(std::string const&);
    void WriteData();

    void SetPrecision(int, int);
    FileName const& DataFilename()    const { return filename_;      }
    Dimension& Dim(Dimension::DimIdxType i) { return Dim_[(int)i];   }
    void SetDim(Dimension::DimIdxType i, Dimension const& d) { Dim_[(int)i]=d; }
    void DataSetNames() const;
    void SetDFLwrite(bool fIn)           { dflWrite_ = fIn;  }
    bool DFLwrite()                const { return dflWrite_; }
  private:
    struct DataFileToken {
      DataFormatType Type;
      const char* Key;
      const char* Description;
      const char* Extension;
      DataIO::AllocatorType Alloc;
    };
    static const DataFileToken DataFileArray[];
    typedef const DataFileToken* TokenPtr;

    static DataIO* AllocDataIO(DataFormatType);
    static DataIO* DetectFormat(std::string const&, DataFormatType&);
    //static DataFormatType DataFormat(std::string const&);

    int debug_;
    int dimension_;         ///< The dimension of all sets in the DataFile.
    DataFormatType dfType_; ///< Format of data in DataFile.
    bool dflWrite_;         ///< If true, write file when DataFileList::WriteAllDF called.
    bool isInverted_;
    bool setDataSetPrecision_; ///< If true set default precision of incoming DataSets.
    int default_width_;
    int default_precision_;
    DataSetList SetList_;      ///< Array of pointers to associated DataSets.
    DataIO* dataio_;           ///< DataIO object for this DataFormatType.
    FileName filename_;        ///< DataFile file name.
    DataIO::DimArray Dim_;     ///< Default file dimensions.
};
#endif
