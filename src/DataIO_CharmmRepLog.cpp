#include <cstdio> //sscanf
#include "DataIO_CharmmRepLog.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"
#include "BufferedLine.h"
#include "DataSet_RemLog.h"

// CONSTRUCTOR
DataIO_CharmmRepLog::DataIO_CharmmRepLog()
{
  SetValid( DataSet::REMLOG );
}

bool DataIO_CharmmRepLog::ID_DataFormat(CpptrajFile& infile) {
  // Assume file set up for read
  if (infile.OpenFile()) return false;
  // Read first two lines
  bool isrepd = false;
  ArgList line1( infile.GetLine() );
  if (line1.Nargs() == 4 &&
      line1[1] == "Replica" &&
      line1[2] == "Exchange")
  {
    const char* ptr = infile.NextLine();
    if (ptr != 0 && ptr[0] == 'R' && ptr[1] == 'E' && ptr[2] == 'X' && ptr[3] == '>')
      isrepd = true;
  }
  infile.CloseFile();
  return isrepd;
}

// DataIO_CharmmRepLog::ReadHelp()
void DataIO_CharmmRepLog::ReadHelp() {
  mprintf("\tnrep <#> : Number of replicas.\n");
}

// DataIO_CharmmRepLog::processReadArgs()
int DataIO_CharmmRepLog::processReadArgs(ArgList& argIn) {
  nrep_ = argIn.getKeyInt("nrep", 0);
  return 0;
}

// DataIO_CharmmRepLog::ReadData()
int DataIO_CharmmRepLog::ReadData(FileName const& fnameIn, 
                            DataSetList& datasetlist, std::string const& dsname)
{
  return ReadReplogArray(fnameIn, datasetlist, dsname);
}

int DataIO_CharmmRepLog::ReadReplogArray(FileName const& fnameIn,
                            DataSetList& datasetlist, std::string const& dsname)
{
  // Expect 1 log for each replica with format <name>_<ext>
  size_t pos = fnameIn.Base().rfind('_');
  std::string prefix = fnameIn.Base().substr(0, pos+1);
  if (debug_ > 0) mprintf("DEBUG: Log file prefix= '%s'\n", prefix.c_str());
  // Search for replica logs from 0 to nrep-1
  typedef std::vector<FileName> Narray;
  Narray Fnames;
  if (nrep_ == 0) {
    mprintf("\tSearching for replica logs with prefix '%s'\n", prefix.c_str());
    FileName fname( fnameIn.DirPrefix() + prefix + integerToString(0) );
    while (File::Exists(fname)) {
      nrep_++;
      Fnames.push_back( fname );
      fname = FileName( fnameIn.DirPrefix() + prefix + integerToString(nrep_) );
    }
  } else {
    for (int i = 0; i != nrep_; i++) {
      FileName fname( fnameIn.DirPrefix() + prefix + integerToString(i) );
      if (!File::Exists(fname)) {
        mprinterr("Error: File '%s' not found.\n", fname.full());
        return 1;
      }
      Fnames.push_back( fname );
    }
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
    /*
    ------------- Replica Exchange ------------
    REX>EXCHANGE =          1  Step =       500
    REX>REPL     =     0  Temp = 313.000  Epot =      -39351.15633548
    REX>NEIGHBOR =    -1  Temp =   0.000  Epot =           0.00000000
    REX>ORIGINAL TAG     0 NEW TAG     0
    REX>PROB     =  0.00000 Rand =  0.52102 Tscale =  1.0000 Success = F
    ------------- Replica Exchange End --------
    */
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
}
