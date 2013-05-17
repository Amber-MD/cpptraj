#ifndef INC_DATAIO_GNUPLOT_H
#define INC_DATAIO_GNUPLOT_H
#include "DataIO.h"
// Class: DataIO_Gnuplot
/// Read/write gnuplot data files.
class DataIO_Gnuplot : public DataIO {
  public:
    DataIO_Gnuplot();
    static DataIO* Alloc() { return (DataIO*)new DataIO_Gnuplot(); }
    int ReadData(std::string const&,DataSetList&) { return 1; }
    int processWriteArgs(ArgList&);
    int WriteData(std::string const&,DataSetList const&,DimArray const&);
    int WriteDataInverted(std::string const&,DataSetList const&,DimArray const&) { return 1; }
    int WriteData2D( std::string const&,DataSet const&, DimArray const&);
    int WriteData3D(std::string const&, DataSet const&, DimArray const&) { return 1; }
    bool ID_DataFormat(CpptrajFile&) { return false; }
  private:
    CpptrajFile file_;
    typedef std::vector<std::string> LabelArray;
    LabelArray Xlabels_;
    LabelArray Ylabels_;
    LabelArray Zlabels_;

    enum PM3D_OPT { OFF = 0, ON, MAP, C2C };
    static const char* BasicPalette[];
    PM3D_OPT pm3d_;
    bool printLabels_; 
    bool useMap_;
    bool jpegout_;
    bool binary_;
    bool writeHeader_;

    static LabelArray LabelArg(std::string const&);
    std::string Pm3d(size_t);
    void WriteRangeAndHeader(Dimension const&, size_t, Dimension const&, size_t,
                             std::string const&);
    void Finish();
    void JpegOut(size_t,size_t);
    void WriteDefinedPalette(int);
    int WriteDataAscii(std::string const&,DataSetList const&, DimArray const&);
    int WriteDataBinary(std::string const&,DataSetList const&, DimArray const&);
};
#endif
