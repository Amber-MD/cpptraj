#include "Exec_Precision.h"
#include "CpptrajStdio.h"

void Exec_Precision::Help() const {
  mprintf("\t{<filename> | <dataset arg>} [<width>] [<precision>]\n"
          "  Set precision for all datasets in datafile <filename> or dataset(s)\n"
          "  specified by <dataset arg> to <width>.<precision>. If width/precision\n"
          "  is not specified then default to 12.4\n");
}

Exec::RetType Exec_Precision::Execute(CpptrajState& State, ArgList& argIn) {
  // Next string is DataSet(s)/DataFile that command pertains to.
  std::string name1 = argIn.GetStringNext();
  if (name1.empty()) {
    mprinterr("Error: No filename/setname given.\n");
    return CpptrajState::ERR;
  }
  // This will break if dataset name starts with a digit...
  int width = argIn.getNextInteger(12);
  if (width < 1) {
    mprintf("Error: Cannot set width < 1 (%i).\n", width);
    return CpptrajState::ERR;
  }
  int precision = argIn.getNextInteger(4);
  if (precision < 0) precision = 0;
  DataFile* df = State.DFL().GetDataFile(name1);
  if (df != 0) {
    mprintf("\tSetting precision for all sets in %s to %i.%i\n", df->DataFilename().base(),
            width, precision);
    df->SetDataFilePrecision(width, precision);
  } else {
    State.DSL().SetPrecisionOfDataSets( name1, width, precision );
  }
  return CpptrajState::OK;
}
