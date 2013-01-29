#ifndef INC_DATAIO_GNUPLOT_H
#define INC_DATAIO_GNUPLOT_H
#include "DataIO.h"
// Class: DataIO_Gnuplot
/// Read/write gnuplot data files.
class DataIO_Gnuplot : public DataIO {
  public:
    DataIO_Gnuplot();

    int processWriteArgs(ArgList &);
    int WriteData(std::string const&,DataSetList&);
    //int WriteDataInverted(DataSetList&);
    int WriteData2D( std::string const&,DataSet& );
  private:
    std::string y_label_;
    double ymin_;
    double ystep_;
    CpptrajFile file_;

    std::vector<std::string> Xlabels_;
    std::vector<std::string> Ylabels_;
    std::vector<std::string> Zlabels_;

    enum PM3D_OPT { OFF = 0, ON, MAP, C2C };
    static const char* BasicPalette[];
    PM3D_OPT pm3d_;
    bool printLabels_; 
    bool useMap_;
    bool jpegout_;
    bool binary_;

    void LabelArg(std::vector<std::string>&, std::string const&);
    std::string Pm3d();
    void WriteRangeAndHeader(double, double, std::string const&);
    void Finish();
    void JpegOut(int,int);
    void WriteDefinedPalette(int);
    int WriteDataAscii(std::string const&,DataSetList&);
    int WriteDataBinary(std::string const&,DataSetList&);
};
#endif
