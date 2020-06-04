#include "DataIO_Peaks.h"
#include "CpptrajStdio.h"
#include "DataSet_Vector_Scalar.h"

/// CONSTRUCTOR
DataIO_Peaks::DataIO_Peaks()
{
  SetValid( DataSet::VECTOR_SCALAR );
}

// DataIO_Peaks::ID_DataFormat()
bool DataIO_Peaks::ID_DataFormat(CpptrajFile& infile)
{

  return false;
}

// DataIO_Peaks::ReadHelp()
void DataIO_Peaks::ReadHelp()
{

}

// DataIO_Peaks::processReadArgs()
int DataIO_Peaks::processReadArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_Peaks::ReadData()
int DataIO_Peaks::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{

  return 1;
}

// DataIO_Peaks::WriteHelp()
void DataIO_Peaks::WriteHelp()
{

}

// DataIO_Peaks::processWriteArgs()
int DataIO_Peaks::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_Peaks::WriteData()
int DataIO_Peaks::WriteData(FileName const& fname, DataSetList const& dsl)
{
  if (dsl.size() > 1)
    mprintf("Warning: Writing multiple sets to peak file may result in invalid format.\n");
  CpptrajFile outfile;
  // Loop over sets. Only write VECTOR_SCALAR sets.
  for (DataSetList::const_iterator it = dsl.begin(); it != dsl.end(); ++it)
  {
    if ((*it)->Type() != DataSet::VECTOR_SCALAR) {
      mprintf("Warning: Set '%s' is not vector with scalar, cannot be used for peaks file.\n",
              (*it)->legend());
    } else {
      DataSet_Vector_Scalar const& ds = static_cast<DataSet_Vector_Scalar const&>( *(*it) );
      if (ds.Size() > 0) {
        // Only open the file when a valid set is found to match old Action_Volmap behavior.
        if (!outfile.IsOpen()) {
          if (outfile.OpenWrite( fname )) {
            mprinterr("Error: Could not open %s for write.\n", fname.full());
            return 1;
          }
        }
        outfile.Printf("%zu\n\n", ds.Size());
        for (unsigned int i = 0; i < ds.Size(); i++)
        {
          Vec3 const& vxyz = ds.Vec(i);
          outfile.Printf("C %16.8f %16.8f %16.8f %16.8f\n",
                         vxyz[0], vxyz[1], vxyz[2], ds.Val(i));
        }
      }
    }
  } // END loop over sets.

  return 0;
}
