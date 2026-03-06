#include "DataIO_AmberFrcmod.h"
#include "CpptrajStdio.h"
#include "DataSet_Parameters.h"
#include "AmberParamFile.h"

/// CONSTRUCTOR
DataIO_AmberFrcmod::DataIO_AmberFrcmod()
{

}

// DataIO_AmberFrcmod::ID_DataFormat()
bool DataIO_AmberFrcmod::ID_DataFormat(CpptrajFile& infile)
{

  return false;
}

// DataIO_AmberFrcmod::ReadHelp()
void DataIO_AmberFrcmod::ReadHelp()
{

}

// DataIO_AmberFrcmod::processReadArgs()
int DataIO_AmberFrcmod::processReadArgs(ArgList& argIn)
{

  return 0;
}

int DataIO_AmberFrcmod::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname) {
  return ReadData(fname, dsl, dsname, 0);
}

// DataIO_AmberFrcmod::ReadData()
int DataIO_AmberFrcmod::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname, DataSet_LeapOpts* leapopts)
{
  ClearAddedByMe();
  // Allocate data set
  MetaData md( dsname );
  DataSet* ds = dsl.CheckForSet( md );
  if (ds != 0) {
    if (ds->Type() != DataSet::PARAMETERS) {
      mprinterr("Error: Set '%s' does not have parameters, cannot append.\n", ds->legend());
      return 1;
    }
    mprintf("\tAdding to existing set %s\n", ds->legend());
  } else {
    ds = dsl.AddSet( DataSet::PARAMETERS, md );
    if (ds == 0) return 1;
  }
  AddedByMe( ds );
  DataSet_Parameters& prm = static_cast<DataSet_Parameters&>( *ds );

  AmberParamFile infile;
  infile.SetAmberParamDebug( debug_ );
  infile.SetDefaults( leapopts );
  if (infile.ReadFrcmod(prm, fname)) {
    mprinterr("Error: Could not read Amber frcmod file '%s'\n", fname.full());
    return 1;
  }

  return 0;
}

// DataIO_AmberFrcmod::WriteHelp()
void DataIO_AmberFrcmod::WriteHelp()
{

}

// DataIO_AmberFrcmod::processWriteArgs()
int DataIO_AmberFrcmod::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_AmberFrcmod::WriteData()
int DataIO_AmberFrcmod::WriteData(FileName const& fname, DataSetList const& dsl)
{

  return 1;
}
