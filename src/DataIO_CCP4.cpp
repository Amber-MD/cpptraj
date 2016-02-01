#include "DataIO_CCP4.h"
#include "CpptrajStdio.h"

// DataIO_CCP4::ID_DataFormat()
bool DataIO_CCP4::ID_DataFormat( CpptrajFile& infile ) {
  bool isCCP4 = false;
  if (!infile.OpenFile()) {
    char MAP[4];
    if (infile.Seek(52 * sizeof(int)) == 0) {
      infile.Read( MAP, 4*sizeof(char) );
      isCCP4 = (MAP[0] == 'M' && MAP[1] == 'A' &&
                MAP[2] == 'P' && MAP[3] == ' ');
    }
    infile.CloseFile();
  }
  return isCCP4;
}

// DataIO_CCP4::ReadData()
int DataIO_CCP4::ReadData(FileName const& fname,
                            DataSetList& datasetlist, std::string const& dsname)
{
  return 1;
}

// DataIO_CCP4::WriteHelp()
void DataIO_CCP4::WriteHelp() {

}

// DataIO_CCP4::processWriteArgs()
int DataIO_CCP4::processWriteArgs(ArgList& argIn) {
  return 0;
}

// DataIO_CCP4::WriteData()
int DataIO_CCP4::WriteData(FileName const& fname, DataSetList const& setList)
{
  return 1;
}
