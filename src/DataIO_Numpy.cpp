#include "DataIO_Numpy.h"
#include "CpptrajStdio.h"

#include "DataSet_Coords.h"

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

/** Read 2d numpy array as COORDS set */
int DataIO_Numpy::read_data_as_coords(std::string const& dsname, DataSetList& dsl,
                                      std::vector<double> const& data,
                                      unsigned long nframes,
                                      unsigned long ncoords)
const
{
  if (ncoords % 3 != 0) {
    mprinterr("Internal Error: DataIO_Numpy::read_data_as_coords(): # coords %lu not divisible by 3.\n", ncoords);
    return 1;
  }
  unsigned long natoms = ncoords / 3;

  DataSet* ds = dsl.AddSet( DataSet::COORDS, MetaData(dsname) );
  if (ds == 0) {
    mprinterr("Error: Could not allocate COORDS set with name '%s'\n", dsname.c_str());
    return 1;
  }
  DataSet_Coords* coords = static_cast<DataSet_Coords*>( ds );

  // Create a pseudo topology
  Topology top;
  for (unsigned long iat = 0; iat != natoms; iat++)
    top.AddTopAtom( Atom("CA","C"), Residue("XXX", 1, ' ', ' ') );
  top.CommonSetup(false, false);
  top.Summary();

  // Set up COORDS set
  if (coords->CoordsSetup( top, CoordinateInfo(Box(), false, false, false) ) ) {
    mprinterr("Error: Could not set up COORDS set '%s'\n", coords->legend());
    return 1;
  }
  if (coords->Allocate(DataSet::SizeArray(1, nframes))) {
    mprinterr("Error: Could not allocate COORDS set '%s'\n", coords->legend());
    return 1;
  }
  Frame frm = coords->AllocateFrame();

  std::vector<double>::const_iterator it = data.begin();
  for (unsigned long ifrm = 0; ifrm != nframes; ifrm++) {
    frm.ClearAtoms();
    for (unsigned long iat = 0; iat != natoms; iat++) {
      double xyz[3];
      xyz[0] = *(it++);
      xyz[1] = *(it++);
      xyz[2] = *(it++);
      frm.AddXYZ( xyz );
    }
    coords->AddFrame( frm );
  }
      

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

  // Convert to coordinates. Expect 2 dimensions, frame and ncoords
  if (dataIn.shape.size() != 2) {
    mprinterr("Error: Shape of numpy array is not 2 (%zu)\n", dataIn.shape.size());
    return 1;
  }
  if (dataIn.fortran_order) {
    mprinterr("Error: Cannot yet process numpy array in fortran order.\n");
    return 1;
  }
  int err = read_data_as_coords(dsname, dsl, dataIn.data, dataIn.shape[0], dataIn.shape[1]);
  if (err != 0) {
    mprinterr("Error: Could not convert numpy array to COORDS.\n");
    return 1;
  }

  return 0;
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
