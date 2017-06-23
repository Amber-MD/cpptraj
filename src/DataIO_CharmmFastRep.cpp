#include <cstdio> //sscanf
#include "DataIO_CharmmFastRep.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"
#include "BufferedLine.h"
#include "DataSet_RemLog.h"

// CONSTRUCTOR
DataIO_CharmmFastRep::DataIO_CharmmFastRep()
{
  SetValid( DataSet::REMLOG );
}

bool DataIO_CharmmFastRep::ID_DataFormat(CpptrajFile& infile) {
  // Assume file set up for read
  if (infile.OpenFile()) return false;
  std::string line = infile.GetLine();
  infile.CloseFile();
  return (line.compare(0,64,"# replica temp. ener. neighbor ntemp nene prob p success? newrep")==0);
}

// DataIO_CharmmFastRep::ReadHelp()
void DataIO_CharmmFastRep::ReadHelp() { }

// DataIO_CharmmFastRep::processReadArgs()
int DataIO_CharmmFastRep::processReadArgs(ArgList& argIn) {
  return 0;
}

inline static int ErrEOF(int line) {
  mprinterr("Error: Unexpected EOF at line %i\n", line);
  return 1;
}

// DataIO_CharmmFastRep::ReadData()
int DataIO_CharmmFastRep::ReadData(FileName const& fnameIn, 
                            DataSetList& datasetlist, std::string const& dsname)
{
  // First determine how many replicas there are.
  BufferedLine infile;
  if (infile.OpenFileRead(fnameIn)) return 1;
  // Skip past header
  const char* ptr = infile.Line();
  // Skip past first Exchange line
  ptr = infile.Line();
  // Position at first replica
  ptr = infile.Line();
  // Sanity check
  if (ptr == 0) return ErrEOF(infile.LineNumber());
  // Get first replica number. Make buffer big enough to avoid overruns.
  char repnum[81];
  int nreps = 0;
  while (ptr != 0 && sscanf(ptr, "%s", repnum) == 1) {
    if (!validInteger(repnum)) break;
    nreps++;
    ptr = infile.Line();
  }
  mprintf("\t%i replicas.\n", nreps);
std::vector<int> CoordinateIndices( nreps );
  // Allocate replica log DataSet
  DataSet* ds = 0;
  if (!dsname.empty()) ds = datasetlist.CheckForSet( dsname );
  if (ds == 0) {
    // New set
    ds = datasetlist.AddSet( DataSet::REMLOG, dsname, "remlog" );
    if (ds == 0) return 1;
    ReplicaDimArray DimTypes;
    DimTypes.AddRemdDimension( ReplicaDimArray::TEMPERATURE );
    ((DataSet_RemLog*)ds)->AllocateReplicas(nreps, DimTypes, 1, false, debug_);
    for (int repidx = 0; repidx != nreps; repidx++)
      CoordinateIndices[repidx] = repidx+1;
  } else {
    if (ds->Type() != DataSet::REMLOG) {
      mprinterr("Error: Set '%s' is not replica log data.\n", ds->legend());
      return 1;
    }
    if ((int)ds->Size() != nreps) {
      mprinterr("Error: Replica log data '%s' is set up for %zu replicas,"
                " current # replicas is %i\n", ds->legend(), ds->Size(),
                nreps);
      return 1;
    }
    mprintf("\tReading final coordinate indices from last frame of existing set.\n");
    for (int repidx = 0; repidx < nreps; repidx++)
      CoordinateIndices[repidx] = ((DataSet_RemLog*)ds)->LastRepFrame(repidx).CoordsIdx();
  }
  mprintf("\tInitial coordinate indices:");
  for (std::vector<int>::const_iterator c = CoordinateIndices.begin();
                                        c != CoordinateIndices.end(); ++c)
    mprintf(" %i", *c);
  mprintf("\n");
  // Close and reopen the file.
  infile.CloseFile();
  if (infile.OpenFileRead(fnameIn)) return 1;
  ptr = infile.Line();
  int nexch = 0;
  bool needs_trim = false;
  DataSet_RemLog& ensemble = static_cast<DataSet_RemLog&>( *ds );
  while (ptr != 0) {
    // Scan to next '# Exchange'
    while (ptr != 0 && !(ptr[0] == '#' && ptr[2] == 'E'))
      ptr = infile.Line();
    if (ptr == 0) break;
    nexch++;
    int crdidx, nbridx;
    double ourtemp, ourpe, nbrtemp, nbrpe;
    char result;
    for (int ourrep = 0; ourrep != nreps; ourrep++) {
      ptr = infile.Line();
      if (sscanf(ptr, "%*2i %12lf %15lf %2i %12lf %15lf %*5f %*5f %c %2i",
                 &ourtemp, &ourpe, &nbridx, &nbrtemp, &nbrpe, &result, &crdidx) != 7)
      {
        needs_trim = true;
        break;
      } 
      // Determine neighbor by difference in temperature
      int nbrrep;
      if (nbridx == -1)
        nbrrep = -1;
      else if (nbrtemp > ourtemp)
        nbrrep = ourrep + 1;
      else
        nbrrep = ourrep - 1; 
//      mprintf("%8i %3i %8.2f %12.4f %3i %8.2f %12.4f %3i %c\n",
//              nexch, ourrep+1, ourtemp, ourpe, nbrrep+1, nbrtemp, nbrpe, crdidx, result);
      ensemble.AddRepFrame( ourrep,
                            DataSet_RemLog::
                            ReplicaFrame( ourrep+1, nbrrep+1, CoordinateIndices[crdidx-1], 0,
                                          (result=='T'), ourtemp, ourpe, 0.0 ) );
    }
  }
  infile.CloseFile();
  mprintf("\tRead %i exchanges.\n", nexch);
  if (needs_trim) ensemble.TrimLastExchange();
  if (debug_ > 1) ensemble.PrintReplicaStats();

  return 0;
}
