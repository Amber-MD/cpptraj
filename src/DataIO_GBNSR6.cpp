#include "DataIO_GBNSR6.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include "AmberEterm.h"
#include "DataSet_double.h"
#include <cstring>

/// CONSTRUCTOR
DataIO_GBNSR6::DataIO_GBNSR6()
{

}

// DataIO_GBNSR6::ID_DataFormat()
bool DataIO_GBNSR6::ID_DataFormat(CpptrajFile& infile)
{

  return false;
}

// DataIO_GBNSR6::ReadHelp()
void DataIO_GBNSR6::ReadHelp()
{

}

// DataIO_GBNSR6::processReadArgs()
int DataIO_GBNSR6::processReadArgs(ArgList& argIn)
{

  return 0;
}

static inline bool validTerm(int i) {
  using namespace Cpptraj;
  AmberEterm::FieldType etype = (AmberEterm::FieldType)i;
  if (etype == AmberEterm::ECAVITY) return true;
  if (etype == AmberEterm::ETOT) return true;
  if (etype == AmberEterm::EEL14) return true;
  if (etype == AmberEterm::EEL) return true;
  if (etype == AmberEterm::EGB) return true;
  if (etype == AmberEterm::ESURF) return true;
  return false;
}
  

// DataIO_GBNSR6::ReadData()
int DataIO_GBNSR6::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  mprintf("\tReading from GBNSR6 output file %s\n", fname.base());
  BufferedLine infile;
  if (infile.OpenFileRead( fname ) != 0) {
    mprinterr("Error: Could not open file '%s'\n", fname.full());
    return 1;
  }
  using namespace Cpptraj;
  AmberEterm AEterm;
  AmberEterm::Darray Energy = AEterm.AllocEnergyArray();
  std::vector<bool> EnergyExists = AEterm.AllocExistsArray();
  DataSetList::DataListType inputSets(AmberEterm::NenergyTerms(), 0);
  const char* ptr = infile.Line();
  enum PhaseType { UNKNOWN=0, INPUT, RESULTS };
  PhaseType Phase = UNKNOWN;
  int frame = 0;
  while (ptr != 0) {
    ArgList argline(ptr);
    if (argline.Nargs() > 0) {
      if (Phase == UNKNOWN) {
        if (argline.Nargs() >= 3 && argline[0] == "Here" && argline[1] == "is" && argline[2] == "the")
          Phase = INPUT;
        else if (argline.Nargs() >= 2 && argline[0] == "3." && argline[1] == "RESULTS")
          Phase = RESULTS;
      } else if (Phase == INPUT) {
        if (strncmp(ptr, "-----", 5) == 0) {
          Phase = UNKNOWN;
        } else {
          if (debug_ > 0) mprintf("DEBUG: [Input] %s\n", ptr);
        }
      } else if (Phase == RESULTS) {
        if (argline[0] == "Maximum")
          Phase = UNKNOWN;
        else if (argline[0] == "Cavity") {
          Energy[AmberEterm::ECAVITY] = argline.getKeyDouble("energy", 0);
          if (debug_ > 0) mprintf("DEBUG: Cavity term: %f\n", Energy[AmberEterm::ECAVITY]);
          EnergyExists[AmberEterm::ECAVITY] = true;
        } else if (strncmp(ptr," ----", 5) == 0) {
          // END frame - store all energies present
          for (int i = 0; i < AmberEterm::NenergyTerms(); i++) {
            if (validTerm(i) && EnergyExists[i]) {
              if (inputSets[i] == 0) {
                MetaData md( dsname, AmberEterm::Ename(i) );
                md.SetLegend( dsname + "_" + AmberEterm::Ename(i) );
                inputSets[i] = new DataSet_double();
                inputSets[i]->SetMeta( md );
              }
              // Since energy terms can appear and vanish over the course of the
              // mdout file, resize if necessary.
              if (frame > (int)inputSets[i]->Size())
                ((DataSet_double*)inputSets[i])->Resize( frame );
              ((DataSet_double*)inputSets[i])->AddElement( Energy[i] );
            }
          }
          frame++;
        } else if (strncmp(ptr, "-----", 5) != 0) {
          if (debug_ > 0) mprintf("DEBUG: [Results] %s\n", ptr);
          if (AEterm.GetAmberEterms(ptr, Energy, EnergyExists))
            mprintf("Warning: Issue parsing line %i\n", infile.LineNumber());
        }
      }
    } // END nargs > 0
    ptr = infile.Line();
  }
  DataSetList::Darray TimeVals(1, 0);
  if (dsl.AddOrAppendSets( "Set", TimeVals, inputSets )) return 1;

  return 0;
}

// DataIO_GBNSR6::WriteHelp()
void DataIO_GBNSR6::WriteHelp()
{

}

// DataIO_GBNSR6::processWriteArgs()
int DataIO_GBNSR6::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_GBNSR6::WriteData()
int DataIO_GBNSR6::WriteData(FileName const& fname, DataSetList const& dsl)
{

  return 1;
}
