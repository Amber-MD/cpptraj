#include "DataIO_AmberEne.h"
#include "BufferedLine.h"
#include "CpptrajStdio.h"
#include <cstdio> // sscanf

/// CONSTRUCTOR
DataIO_AmberEne::DataIO_AmberEne()
{

}

// DataIO_AmberEne::ID_DataFormat()
bool DataIO_AmberEne::ID_DataFormat(CpptrajFile& infile)
{
  if (infile.OpenFile()) return false;
  bool isAmberEne = false;
  std::string line = infile.GetLine();
  ArgList lineArgs( line, " " );
  if (lineArgs.Nargs() > 2) {
    if (lineArgs[0] == "L0" && lineArgs[1] == "Nsteps")
      isAmberEne = true;
  }
  return isAmberEne;
}

// DataIO_AmberEne::ReadHelp()
void DataIO_AmberEne::ReadHelp()
{

}

// DataIO_AmberEne::processReadArgs()
int DataIO_AmberEne::processReadArgs(ArgList& argIn)
{

  return 0;
}

// ----------------------------------------------- 
/// Class to hold characteristics of an Amber energy file line
class AmberEneLine {
  public:
    /// CONSTRUCTOR
    AmberEneLine() {}
    /// Set up from given line
    int Setup(const char* ptr, unsigned int lineNumber, DataSetList& dsl,
              std::string const& dsname, DataSet::DataType fpType)
    {
      ArgList headerArgs( ptr, " \n\r" );
      if (headerArgs.Nargs() < 1) {
        mprinterr("Error: No columns detected at line %u of Amber energy file.\n",
                  lineNumber);
        return 1;
      }
      //fmt_.assign("%*s");
      // Output sets and sscanf read format
      for (int iarg = 1; iarg < headerArgs.Nargs(); iarg++) {
      // Determine expected type
        DataSet::DataType setType;
        if (headerArgs[iarg] == "Nsteps")
          setType = DataSet::INTEGER;
        else
          setType = fpType;
        MetaData md(dsname, headerArgs[iarg]);
        DataSet* ds = dsl.CheckForSet( md );
        if (ds == 0) {
          // Creating new set
          ds = dsl.AddSet( setType, md );
          if (ds == 0) {
            mprinterr("Error: Could not allocate set '%s[%s]'\n",
                      dsname.c_str(), headerArgs[iarg].c_str());
            return 1;
          }
        } else {
          // Append to existing. Check type.
          if (setType != ds->Type()) {
            mprinterr("Error: Cannot append to set '%s'. Type is '%s', expected '%s'\n",
                      ds->legend(), DataSet::description(ds->Type()),
                      DataSet::description(setType));
            return 1;
          }
        }
        outsets_.push_back( ds );
        //if (setType == DataSet::INTEGER)
        //  fmt_.append(" %i");
        //else
        //  fmt_.append(" %lf");
      } // END loop over header labels in this line
      return 0;
    } // END Setup()
    /// Print details to stdout
    void Print() const {
      //mprintf("Format: '%s'\n", fmt_.c_str());
      mprintf("\tSets:");
      for (std::vector<DataSet*>::const_iterator it = outsets_.begin();
                                                 it != outsets_.end(); ++it)
        mprintf(" %s", (*it)->legend());
      mprintf("\n");
    }
    /// \return # of values this line contains
    unsigned int Nvals() const { return outsets_.size(); }
    /// \return iterator to beginning of sets
    std::vector<DataSet*>::const_iterator begin() const { return outsets_.begin(); }
    /// \return iterator to end of sets
    std::vector<DataSet*>::const_iterator end()   const { return outsets_.end(); }
  private:
    //std::string fmt_;               ///< scanf format for reading this line
    std::vector<DataSet*> outsets_; ///< Output sets for this line
};
// ----------------------------------------------- 

// DataIO_AmberEne::ReadData()
int DataIO_AmberEne::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  BufferedLine infile;
  if (infile.OpenFileRead(fname)) return 1;
  const char* ptr = infile.Line();
  // Should be at first line of the header
  //mprintf("DEBUG: %i '%s'\n", infile.LineNumber(), ptr);
  if (ptr[0] != 'L' || ptr[1] != '0') {
    mprinterr("Error: 'L0' not found in Amber energy file '%s'\n", fname.full());
    return 1;
  }

  typedef std::vector<AmberEneLine> LineArray;
  LineArray Lines;

  // Read the header
  unsigned int maxVals = 0;
  const DataSet::DataType fpType = DataSet::DOUBLE;
  while (ptr != 0) {
   if (infile.LineNumber() > 1 && ptr[0] == 'L' && ptr[1] == '0') {
      // At first line of data
      break;
    }
    Lines.push_back( AmberEneLine() );
    if (Lines.back().Setup(ptr, infile.LineNumber(), dsl, dsname, fpType))
      return 1;
    Lines.back().Print();
    maxVals = std::max( maxVals, Lines.back().Nvals() );

    ptr = infile.Line();
  }

  // Sanity check
  if (Lines.empty()) {
    mprinterr("Error: No headers read.\n");
    return 1;
  }

  // Read data
  int ndata = 0;
  while (ptr != 0) {
    for (LineArray::const_iterator it = Lines.begin(); it != Lines.end(); ++it)
    {
      int nvals = infile.TokenizeLine(" \n\r");
      // Subtract the LX column
      nvals--;
      if ( (nvals < 1) || ((unsigned int)nvals != it->Nvals()) ) {
        mprinterr("Error: Line %i expected %u values, got %i\n",
                  infile.LineNumber(), it->Nvals(), nvals);
        it->Print();
        return 1;
      }
      infile.NextToken();
      // Loop over remaining columns
      for (std::vector<DataSet*>::const_iterator ds = it->begin();
                                                 ds != it->end(); ++ds)
      {
        if ((*ds)->Type() == DataSet::INTEGER) {
          int ival = atoi( infile.NextToken() );
          (*ds)->Add( ndata, &ival );
        } else {
          double dval = atof( infile.NextToken() );
          (*ds)->Add( ndata, &dval );
        }
      } // END loop over remaining columns
      ptr = infile.Line();
    } // END loop over lines

    ndata++;
  }
  return 0;
}

// DataIO_AmberEne::WriteHelp()
void DataIO_AmberEne::WriteHelp()
{

}

// DataIO_AmberEne::processWriteArgs()
int DataIO_AmberEne::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_AmberEne::WriteData()
int DataIO_AmberEne::WriteData(FileName const& fname, DataSetList const& dsl)
{

  return 1;
}
