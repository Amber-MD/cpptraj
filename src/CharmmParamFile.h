#ifndef INC_CHARMMPARAMFILE_H
#define INC_CHARMMPARAMFILE_H
#include <string>
class ParameterSet;
class BufferedLine;
class FileName;
/// Used to read in CHARMM parameters from CHARMM parameter file.
class CharmmParamFile {
  public:
    CharmmParamFile() {}
    int ReadParams(ParameterSet&, FileName const&, int) const;
    int WriteParams(ParameterSet&, FileName const&, int) const;
  private:
    int ReadInput(std::string&, BufferedLine&) const;
};
#endif
