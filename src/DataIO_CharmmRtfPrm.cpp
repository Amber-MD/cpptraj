#include <cctype> // isspace
#include "DataIO_CharmmRtfPrm.h"
#include "CpptrajStdio.h"
#include "CharmmParamFile.h"
#include "DataSet_Parameters.h"

/// CONSTRUCTOR
DataIO_CharmmRtfPrm::DataIO_CharmmRtfPrm()
{

}

// DataIO_CharmmRtfPrm::ID_DataFormat()
bool DataIO_CharmmRtfPrm::ID_DataFormat(CpptrajFile& infile)
{

  return false;
}

// DataIO_CharmmRtfPrm::ReadHelp()
void DataIO_CharmmRtfPrm::ReadHelp()
{

}

// DataIO_CharmmRtfPrm::processReadArgs()
int DataIO_CharmmRtfPrm::processReadArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_CharmmRtfPrm::ReadData()
int DataIO_CharmmRtfPrm::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  mprintf("Warning: Currently only CHARMM parameters will be read from this file.\n");
 
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
  DataSet_Parameters& prm = static_cast<DataSet_Parameters&>( *ds );

  CharmmParamFile infile;
  return infile.ReadParams( prm, fname, debug_ );
}

// DataIO_CharmmRtfPrm::WriteHelp()
void DataIO_CharmmRtfPrm::WriteHelp()
{

}

// DataIO_CharmmRtfPrm::processWriteArgs()
int DataIO_CharmmRtfPrm::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_CharmmRtfPrm::WriteData()
int DataIO_CharmmRtfPrm::WriteData(FileName const& fname, DataSetList const& dsl)
{

  return 1;
}
