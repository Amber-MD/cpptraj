#include <cctype> // isspace
#include "DataIO_CharmmRtfPrm.h"
#include "CpptrajStdio.h"
#include "CharmmParamFile.h"
#include "DataSet_Parameters.h"
#include "DataSet_Topology.h"

/// CONSTRUCTOR
DataIO_CharmmRtfPrm::DataIO_CharmmRtfPrm()
{
  SetValid( DataSet::PARAMETERS );
  SetValid( DataSet::TOPOLOGY );
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
  ParameterSet pout;
  ParameterSet::UpdateCount ucount;
  for (DataSetList::const_iterator it = dsl.begin(); it != dsl.end(); ++it)
  {
    if ( (*it)->Type() == DataSet::PARAMETERS )
    {
      mprintf("\tUsing parameter set '%s'\n", (*it)->legend());
      DataSet_Parameters const& dsPrm = static_cast<DataSet_Parameters const&>( *(*it) );
      pout.UpdateParamSet( dsPrm, ucount, debug_ );
    }
    else if ( (*it)->Type() == DataSet::TOPOLOGY )
    {
      mprintf("\tUsing parameters from topology '%s'\n", (*it)->legend());
      // Convert topology to parameter set
      DataSet_Topology const& dsTop = static_cast<DataSet_Topology const&>( *(*it) );
      pout.UpdateParamSet( dsTop.Top().GetParameters(), ucount, debug_ );
    } else {
      mprintf("Warning: '%s' is not a valid parameter/topology set, skipping.\n",
              (*it)->legend());
    }
  }
  // TODO check if any parameters were added
  CharmmParamFile outfile;
  return outfile.WriteParams( pout, fname, debug_ );
}
