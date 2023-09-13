#include "DataIO_Numpy.h"
#include "CpptrajStdio.h"

#include "libnpy/npy.hpp"

/// CONSTRUCTOR
DataIO_Numpy::DataIO_Numpy()
{

}

// DataIO_Numpy::ID_DataFormat()
bool DataIO_Numpy::ID_DataFormat(CpptrajFile& infile)
{
  if (infile.OpenFile()) return false;
  // Read first 6 bytes
  unsigned char magic[6];
  magic[0] = 0;
  magic[1] = 0;
  magic[2] = 0;
  magic[3] = 0;
  magic[4] = 0;
  magic[5] = 0;
  infile.Read(magic, 6);
  infile.CloseFile();
  if (magic[0] == 93  && magic[1] == 'N' && magic[2] == 'U' &&
      magic[3] == 'M' && magic[4] == 'P' && magic[5] == 'Y')
    return true;
  return false;
}

// DataIO_Numpy::ReadHelp()
void DataIO_Numpy::ReadHelp()
{

}

// DataIO_Numpy::processReadArgs()
int DataIO_Numpy::processReadArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_Numpy::ReadData()
int DataIO_Numpy::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  mprintf("\tReading numpy array from '%s'\n", fname.full());
  npy::npy_data<double> dataIn = npy::read_npy<double>( fname.Full() );

  //std::vector<double> data = dataIn.data;
  //std::vector<unsigned long> shape = dataIn.shape;
  //bool fortran_order = dataIn.fortran_order;

  mprintf("\tData array size = %zu\n", dataIn.data.size());
  mprintf("\tShape array size = %zu\n", dataIn.shape.size());
  for (std::vector<unsigned long>::const_iterator it = dataIn.shape.begin();
                                                  it != dataIn.shape.end(); ++it)
    mprintf("\t\t%lu\n", *it);
  mprintf("\tFortran order = %i\n", (int)dataIn.fortran_order);

  return 1;
}

// DataIO_Numpy::WriteHelp()
void DataIO_Numpy::WriteHelp()
{

}

// DataIO_Numpy::processWriteArgs()
int DataIO_Numpy::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_Numpy::WriteData()
int DataIO_Numpy::WriteData(FileName const& fname, DataSetList const& dsl)
{
  mprintf("\tWriting to numpy array '%s'\n", fname.full());
  return 1;
}
