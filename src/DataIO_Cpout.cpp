#include <cstdio>
#include "DataIO_Cpout.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include "DataSet_PH.h"

/// CONSTRUCTOR
DataIO_Cpout::DataIO_Cpout() :
  type_(NONE)
{
  SetValid( DataSet::PH );
}

const char* DataIO_Cpout::FMT_REDOX_ = "Redox potential: %f V";

const char* DataIO_Cpout::FMT_PH_ = "Solvent pH: %f";

// DataIO_Cpout::ID_DataFormat()
bool DataIO_Cpout::ID_DataFormat(CpptrajFile& infile)
{
  bool iscpout = false;
  type_ = NONE;
  if (!infile.OpenFile()) {
    const char* ptr = infile.NextLine();
    if (ptr != 0) {
      float orig_ph;
      if (sscanf(ptr, FMT_REDOX_, &orig_ph) == 1) {
        type_ = REDOX;
      } else if (sscanf(ptr, FMT_PH_, &orig_ph) == 1) {
        type_ = PH;
      }
      if (type_ != NONE) {
        ptr = infile.NextLine();
        int step_size;
        if (ptr != 0) {
          iscpout = (sscanf(ptr, "Monte Carlo step size: %d", &step_size) == 1);
        }
      }
    }
    infile.CloseFile();
  }
  return iscpout;
}

// DataIO_Cpout::ReadHelp()
void DataIO_Cpout::ReadHelp()
{

}

// DataIO_Cpout::processReadArgs()
int DataIO_Cpout::processReadArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_Cpout::ReadData()
int DataIO_Cpout::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  BufferedLine infile;
  if (infile.OpenFileRead( fname )) return 1;
  const char* ptr = infile.Line();

  float orig_ph, time, pHval = 0.0;
  int step, res, state;

  // Determine type if necessary.
  if (type_ == NONE) {
    if (sscanf(ptr, FMT_REDOX_, &orig_ph) == 1) {
      type_ = REDOX;
    } else if (sscanf(ptr, FMT_PH_, &orig_ph) == 1) {
      type_ = PH;
    } else {
      mprinterr("Error: Could not determine CPOUT file type.\n");
      return 1;
    }
    infile.CloseFile();
    infile.OpenFileRead( fname );
  }

  const char* fmt = 0;
  const char* rFmt = 0;
  if (type_ == PH) {
    mprintf("\tConstant pH output file.\n");
    fmt = FMT_PH_;
    rFmt = "Residue %d State: %d pH: %f";
  } else if (type_ == REDOX) {
    mprintf("\tRedOx output file.\n");
    fmt = FMT_REDOX_;
    rFmt = "Residue %d State: %d E: %f V";
  }

  // Allocate ph DataSet
  DataSet* ds = 0;
  if (!dsname.empty()) ds = dsl.CheckForSet(dsname);
  if (ds == 0) {
    // New set
    ds = dsl.AddSet( DataSet::PH, dsname, "ph" );
    if (ds == 0) return 1;
  } else {
    if (ds->Type() != DataSet::PH) {
      mprinterr("Error: Set '%s' is not ph data.\n", ds->legend());
      return 1;
    }
    // TODO check # residues etc
  }
  DataSet_PH* phdata = (DataSet_PH*)ds;

  while (ptr != 0) {
    if (sscanf(ptr, fmt, &orig_ph) == 1) {
      // Full record
      mprintf("DEBUG: pH= %f\n", orig_ph);
      ptr = infile.Line(); // Monte Carlo step size
      ptr = infile.Line(); // Current MD time step
      if (sscanf(ptr,"Time step: %d", &step) != 1) {
        mprinterr("Error: Could not get step.\n");
        return 1;
      }
      mprintf("DEBUG: step= %i\n", step);
      ptr = infile.Line(); // Current time (ps)
      if (sscanf(ptr,"Time: %f", &time) != 1) {
        mprinterr("Error: Could not get time.\n");
        return 1;
      }
      mprintf("DEBUG: time= %f\n", time);
      ptr = infile.Line(); // Residue
    } 
    // delta record or full record header read
    while (sscanf(ptr, rFmt, &res, &state, &pHval) >= 2) {
      mprintf("DEBUG: res= %i state= %i pH= %f\n", res, state, pHval);
      phdata->AddState(res, state);
      ptr = infile.Line();
    }
    ptr = infile.Line();
  }
  infile.CloseFile();

  return 0;
}

// DataIO_Cpout::WriteHelp()
void DataIO_Cpout::WriteHelp()
{

}

// DataIO_Cpout::processWriteArgs()
int DataIO_Cpout::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_Cpout::WriteData()
int DataIO_Cpout::WriteData(FileName const& fname, DataSetList const& dsl)
{

  return 1;
}
