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
  mprintf("DEBUG: '%s'\n", line.c_str());
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
  mprintf("DEBUG: %i replicas.\n", nreps);
  // Close and reopen the file.
  infile.CloseFile();
  if (infile.OpenFileRead(fnameIn)) return 1;
  ptr = infile.Line();
  int nexch = 0;
  while (ptr != 0) {
    // Scan to next '# Exchange'
    while (ptr != 0 && !(ptr[0] == '#' && ptr[2] == 'E'))
      ptr = infile.Line();
    if (ptr == 0) break;
    mprintf("%s\n", ptr);
    nexch++;
    int crdidx, nbridx;
    double ourtemp, ourpe, nbrtemp, nbrpe;
    char result;
    for (int ourrep = 0; ourrep != nreps; ourrep++) {
      ptr = infile.Line();
      if (sscanf(ptr, "%2i %12lf %15lf %2i %12lf %15lf %*5f %*5f %c",
                 &crdidx, &ourtemp, &ourpe, &nbridx, &nbrtemp, &nbrpe, &result) != 7)
        return 1;
      // Determine neighbor by difference in temperature
      int nbrrep;
      if (nbridx == -1)
        nbrrep = -1;
      else if (nbrtemp > ourtemp)
        nbrrep = ourrep + 1;
      else
        nbrrep = ourrep - 1; 
      mprintf("%8i %3i %8.2f %12.4f %3i %8.2f %12.4f %3i %c\n",
              nexch, ourrep+1, ourtemp, ourpe, nbrrep+1, nbrtemp, nbrpe, crdidx, result);

      //ensemble.AddRepFrame( ourrep,
      //                      DataSet_RemLog::
      //                      ReplicaFrame( ourrep, nbrrep, CoordinateIndices[crdidx], 0,
      //                                    result, ourtemp, ourpe, 0.0 ) );
    }
  }

      
/*


  // Expect 1 log for each replica with format <name>_<ext>
  size_t pos = fnameIn.Base().rfind('_');
  std::string prefix = fnameIn.Base().substr(0, pos+1);
  if (debug_ > 0) mprintf("DEBUG: Log file prefix= '%s'\n", prefix.c_str());
  // Search for replica logs from 0 to nrep-1
  typedef std::vector<FileName> Narray;
  Narray Fnames;
  for (int i = 0; i != nrep_; i++) {
    FileName fname( fnameIn.DirPrefix() + prefix + integerToString(i) );
    if (!File::Exists(fname)) {
      mprinterr("Error: File '%s' not found.\n", fname.full());
      return 1;
    }
    Fnames.push_back( fname );
  }
  mprintf("\t%zu replica logs.\n", Fnames.size());
  std::vector<int> CoordinateIndices( nrep_ );
  // Allocate replica log DataSet
  DataSet* ds = 0;
  if (!dsname.empty()) ds = datasetlist.CheckForSet( dsname );
  if (ds == 0) {
    // New set
    ds = datasetlist.AddSet( DataSet::REMLOG, dsname, "remlog" );
    if (ds == 0) return 1;
    // FIXME assume temperature for now
    ReplicaDimArray DimTypes;
    DimTypes.AddRemdDimension( ReplicaDimArray::TEMPERATURE );
    ((DataSet_RemLog*)ds)->AllocateReplicas(nrep_, DimTypes, 0, false, debug_);
    for (int repidx = 0; repidx != nrep_; repidx++)
      CoordinateIndices[repidx] = repidx;
  } else {
    if (ds->Type() != DataSet::REMLOG) {
      mprinterr("Error: Set '%s' is not replica log data.\n", ds->legend());
      return 1;
    }
    if ((int)ds->Size() != nrep_) {
      mprinterr("Error: Replica log data '%s' is set up for %zu replicas,"
                " current # replicas is %i\n", ds->legend(), ds->Size(),
                nrep_);
      return 1;
    }
    mprintf("\tReading final coordinate indices from last frame of existing set.\n");
    for (int repidx = 0; repidx < nrep_; repidx++)
      CoordinateIndices[repidx] = ((DataSet_RemLog*)ds)->LastRepFrame(repidx).CoordsIdx();
  }
  mprintf("\tInitial coordinate indices:");
  for (std::vector<int>::const_iterator c = CoordinateIndices.begin();
                                        c != CoordinateIndices.end(); ++c)
    mprintf(" %i", *c);
  mprintf("\n");
  // Loop over replica logs
  DataSet_RemLog& ensemble = static_cast<DataSet_RemLog&>( *ds );
  int total_exchanges = -1;
  bool needs_trim = false;
  bool isHREMD = false;
  for (int i = 0; i != nrep_; i++) {
    mprintf("\t\t%s\n", Fnames[i].full());
    bool warnsgld = false;
    BufferedLine infile;
    if (infile.OpenFileRead( Fnames[i] )) return 1;
    *
    ------------- Replica Exchange ------------
    REX>EXCHANGE =          1  Step =       500
    REX>REPL     =     0  Temp = 313.000  Epot =      -39351.15633548
    REX>NEIGHBOR =    -1  Temp =   0.000  Epot =           0.00000000
    REX>ORIGINAL TAG     0 NEW TAG     0
    REX>PROB     =  0.00000 Rand =  0.52102 Tscale =  1.0000 Success = F
    ------------- Replica Exchange End --------
    *
    const char* ptr = infile.Line();
    int nexch = 0;
    while (ptr != 0) {
      // First line is 'Replica Exchange'
      // Next line is exchange and step
      ptr = infile.Line();
      // Next line is replica, temperature, my PE
      ptr = infile.Line();
      int ourrep;
      double ourtemp;
      double ourpe;
      sscanf(ptr, "%*15c%5i%*9c%7lf%*9c%20lf", &ourrep, &ourtemp, &ourpe);
      if (ourrep < 0 || ourrep >= nrep_) {
        mprinterr("Error: Replica number %i is out of range.\n", ourrep);
        return 1;
      }
      // Next line is neighbor, neighbor temp, neighbor PE
      ptr = infile.Line();
      int nbrrep;
      double nbrtemp;
      double nbrpe;
      sscanf(ptr, "%*15c%5i%*9c%7lf%*9c%20lf", &nbrrep, &nbrtemp, &nbrpe);
      // Neighbor -1 is no neighbor
      if (nbrrep < -1 || nbrrep >= nrep_) {
        mprinterr("Error: Neighbor replica number %i is out of range.\n", nbrrep);
        return 1;
      }
      // Next line may be tag, i.e. coordinate index, but check some other stuff first
      ptr = infile.Line();
      if (ptr[0]=='R' && ptr[1]=='X') {
        // Skip RXSGLD stuff
        if (!warnsgld) {
          mprintf("Warning: Log has RXSGLD info - this info is being skipped.\n");
          warnsgld = true;
        }
        ptr = infile.Line();
      }
      // Check for HREMD info
      double pe_x2 = 0.0;
      if (ptr[0]=='T' && ptr[1]=='H') {
        sscanf(ptr, "%*6c%*12f%*12f%*12f%12lf", &pe_x2);
        ptr = infile.Line();
        if (!isHREMD)
          mprintf("Info: Hamiltonian (THAM) info detected in log.\n");
        isHREMD = true;
      }
      // Now actually at the coordinate index
      int crdidx;
      sscanf(ptr, "%*31c%5i", &crdidx);
      // Next line has result
      ptr = infile.Line();
      bool result = (ptr[67]=='T');
//      mprintf("%8i %3i %6.2f %12.4f %3i %6.2f %12.4f %12.4f %c\n",
//              nexch, ourrep, ourtemp, ourpe, nbrrep, nbrtemp, nbrpe, pe_x2, ptr[67]);
      // FIXME: Need +1 for ourrep and nbrrep?
      ensemble.AddRepFrame( ourrep,
                            DataSet_RemLog::
                            ReplicaFrame( ourrep, nbrrep, CoordinateIndices[crdidx], 0,
                                          result, ourtemp, ourpe, pe_x2 ) );
      // Scan to Replica Exchange End
      while (ptr != 0 && *ptr != '-')
        ptr = infile.Line();
      nexch++;
      // Next line is beginning of next exchange
      ptr = infile.Line();
    }
    infile.CloseFile();
    if (total_exchanges < 0)
      total_exchanges = nexch;
    else if (nexch != total_exchanges) {
      mprintf("Warning: Number of exchanges %i != # in first file %i\n", nexch, total_exchanges);
      needs_trim = true;
    }
  } // End loop over replica logs
  if (isHREMD) ensemble.SetupDimTypes().ChangeRemdDim(0, ReplicaDimArray::HAMILTONIAN);
  if (needs_trim) ensemble.TrimLastExchange();
  if (debug_ > 1) ensemble.PrintReplicaStats();

  return 0;
*/
  return 1;
}
