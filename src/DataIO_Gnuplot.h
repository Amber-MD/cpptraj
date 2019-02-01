#ifndef INC_DATAIO_GNUPLOT_H
#define INC_DATAIO_GNUPLOT_H
#include "DataIO.h"
#include "BufferedLine.h"
/// Read/write gnuplot data files.
class DataIO_Gnuplot : public DataIO {
  public:
    DataIO_Gnuplot();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_Gnuplot(); }
    static void WriteHelp();
    int processReadArgs(ArgList&) { return 0; }
    int ReadData(FileName const&,DataSetList&,std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&,DataSetList const&);
    bool ID_DataFormat(CpptrajFile&);
  private:
    CpptrajFile file_;
    FileName data_fname_; ///< Data file name
    std::string title_;   ///< Plot title (output)
    typedef std::vector<std::string> LabelArray;
    LabelArray Xlabels_;
    LabelArray Ylabels_;
    LabelArray Zlabels_;

    enum PM3D_OPT { OFF = 0, ON, MAP, C2C };
    static const char* BasicPalette[];
    PM3D_OPT pm3d_;
    std::string palette_;
    bool printLabels_; 
    bool useMap_;
    bool jpegout_;
    bool binary_;
    bool writeHeader_;

    int ReadBinaryData(FileName const&,DataSetList&,std::string const&,
                       std::string const&, std::string const&);
    int ReadAsciiHeader(FileName const&,DataSetList&,std::string const&);
    int ReadAsciiData(BufferedLine&, DataSetList&, std::string const&,
                      std::string const&, std::string const&);

    static LabelArray LabelArg(std::string const&);
    int WriteSet2D( DataSet const& );
    std::string Pm3d(size_t);
    void WriteRangeAndHeader(Dimension const&, size_t, Dimension const&, size_t,
                             std::string const&);
    void WriteLabels(LabelArray const&, Dimension const&, char);
    void Finish();
    void JpegOut(size_t,size_t);
    void WriteDefinedPalette(int);
    int WriteSets1D(DataSetList const&);
};
#endif
