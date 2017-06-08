#include <cstdio> //sscanf
#include "DataIO_CharmmRepLog.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"

// CONSTRUCTOR
DataIO_CharmmRepLog::DataIO_CharmmRepLog()
{
  SetValid( DataSet::REMLOG );
}

// NOTE: Must match LogType
const char* DataIO_CharmmRepLog::LogDescription[] = {
  "Unknown", "Temperature", "Hamiltonian", "MultipleDim", "RXSGLD", "pH"
};

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
  mprintf("\tnosearch            : Do not automatically search for MREMD dimension logs.\n"
          "\tdimfile <file>      : remd.dim file for processing MREMD logs.\n"
          "\tcrdidx <crd indices>: Use comma-separated list of indices as the initial\n"
          "\t                      coordinate indices.\n"
          "\tMultiple REM logs may be specified.\n");
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
  int err = 0;
  if (nrep_ > 0)
    err = ReadReplogArray(fnameIn, datasetlist, dsname);
 
  for (Farray::iterator it = Logs_.begin(); it != Logs_.end(); ++it)
    it->CloseFile();
 
  return 1;
}

int DataIO_CharmmRepLog::ReadReplogArray(FileName const& fnameIn,
                            DataSetList& datasetlist, std::string const& dsname)
{
  // Expect 1 log for each replica with format <name>_<ext>
  size_t pos = fnameIn.Base().rfind('_');
  std::string prefix = fnameIn.Base().substr(0, pos+1);
  mprintf("DEBUG: prefix= '%s'\n", prefix.c_str());
  // Search for replica logs from 0 to nrep-1
  for (int i = 0; i != nrep_; i++) {
    FileName fname( prefix + integerToString(i) );
    if (!File::Exists(fname)) {
      mprinterr("Error: File '%s' not found.\n", fname.full());
      return 1;
    }
    mprintf("DEBUG:\t\t%s\n", fname.full());
    bool warnsgld = false;
    BufferedLine infile;
    if (infile.OpenFileRead( fname )) return 1;
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
      // Next line is neighbor, neighbor temp, neighbor PE
      ptr = infile.Line();
      int nbrrep;
      double nbrtemp;
      double nbrpe;
      sscanf(ptr, "%*15c%5i%*9c%7lf%*9c%20lf", &nbrrep, &nbrtemp, &nbrpe);
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
      double pe_x2 = 0.0;
      if (ptr[0]=='T' && ptr[1]=='H') {
        // Log has HREMD info
        sscanf(ptr, "%*6c%*12f%*12f%*12f%12lf", &pe_x2);
        ptr = infile.Line();
      }
      // Now actually at the coordinate index
      int crdidx;
      sscanf(ptr, "%*17c%5i", &crdidx);
      // Next line has result
      ptr = infile.Line();
      bool result = (ptr[67]=='T');
      mprintf("%8i %3i %6.2f %12.4f %3i %6.2f %12.4f %12.4f %c\n",
              nexch++, ourrep, ourtemp, ourpe, nbrrep, nbrtemp, nbrpe, pe_x2, ptr[67]);
      // Next line is Replica Exchange End
      ptr = infile.Line();
      // Next line is beginning of next exchange
      ptr = infile.Line();
    }
    infile.CloseFile();
  }

  return 1;
}
