#include "DataIO_Peaks.h"
#include "CpptrajStdio.h"
#include "DataSet_Vector_Scalar.h"
#include "StringRoutines.h"
#include <cstdio> // sscanf

/// CONSTRUCTOR
DataIO_Peaks::DataIO_Peaks()
{
  SetValid( DataSet::VECTOR_SCALAR );
}

// DataIO_Peaks::ID_DataFormat()
bool DataIO_Peaks::ID_DataFormat(CpptrajFile& infile)
{
  if (infile.OpenFile()) return false;
  std::string firstLine = NoTrailingWhitespace(infile.GetLine());
  std::string secondLine = NoTrailingWhitespace(infile.GetLine());
  std::string thirdLine = infile.GetLine();
  //mprintf("DEBUG: ID peaks.\n");
  //mprintf("'%s'\n", firstLine.c_str());
  //mprintf("'%s'\n", secondLine.c_str());
  //mprintf("'%s'\n", thirdLine.c_str());
  infile.CloseFile();
  if (validInteger(firstLine)) {
    //mprintf("DEBUG: Line 1 int.\n");
    if (secondLine.empty()) {
      //mprintf("DEBUG: Line 2 empty.\n");
      ArgList line3(thirdLine);
      if (line3.Nargs() == 5) {
        if (line3[0] == "C") {
          for (int col = 1; col < 5; col++)
            if (!validDouble(line3[col])) return false;
          //mprintf("DEBUG: Line 3 C and 4 doubles.\n");
          return true;
        }
      }
    }
  }
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
  // Figure out the data set
  unsigned int peakidx = 0;
  DataSet* ds = dsl.CheckForSet( dsname );
  if (ds == 0) {
    // New set
    ds = dsl.AddSet( DataSet::VECTOR_SCALAR, dsname, "peaks" );
    if (ds == 0) return 1;
  } else {
    // Appending.
    if (ds->Type() != DataSet::VECTOR_SCALAR) {
      mprinterr("Error: Set '%s' is not vector with scalar, cannot append peaks to it.\n",
                ds->legend());
      return 1;
    }
    peakidx = ds->Size();
    mprintf("\tAppending peaks to set '%s', offset %u\n", ds->legend(), peakidx);
  }
  // Parse through the peaks file and extract the peaks
  CpptrajFile peakfile;
  if (peakfile.OpenRead(fname)) {
    mprinterr("Error: Could not open %s for reading!\n", fname.full());
    return 1;
  }
  // FORMAT:
  //   Line 1   : # peaks
  int npeaks_in_file = 0;
  const char* ptr = peakfile.NextLine();
  if (ptr == 0) {
    mprinterr("Error: Unexpected EOF when reading # peaks..\n");
    return 1;
  }
  if (sscanf(ptr, "%i", &npeaks_in_file) != 1) {
    mprinterr("Error: Could not read number of peaks in file.\n");
    return 1;
  }
  mprintf("\tAttempting to read %i peaks.\n", npeaks_in_file);
  //   Line 2   : Blank
  ptr = peakfile.NextLine();
  //   Line 3-X : C <X> <Y> <Z> <Val>
  int npeaks_read = 0;
  ptr = peakfile.NextLine();
  double pbuf[4];
  while (ptr != 0) {
    if (sscanf(ptr, "C %lg %lg %lg %lg", pbuf, pbuf+1, pbuf+2, pbuf+3) != 4) {
      mprintf("Warning: Unexpected format for line, skipping: '%s'\n",
              NoTrailingWhitespace(std::string(ptr)).c_str());
    } else {
      npeaks_read++;
      ds->Add(peakidx++, pbuf);
    }
    ptr = peakfile.NextLine();
  }
  peakfile.CloseFile();
  mprintf("\tRead %i peaks.\n", npeaks_read);
  if (npeaks_read != npeaks_in_file)
    mprintf("Warning: Expected %i peaks, actually got %i peaks.\n", npeaks_in_file, npeaks_read);

  return 0;
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
