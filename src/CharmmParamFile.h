#ifndef INC_CHARMMPARAMFILE_H
#define INC_CHARMMPARAMFILE_H
#include <string>
class BufferedLine;
class FileName;
namespace Cpptraj {
namespace Parm {
class ParameterSet;
}
}
/// Used to read in CHARMM parameters from CHARMM parameter file.
class CharmmParamFile {
  public:
    CharmmParamFile() {}
    int ReadParams(Cpptraj::Parm::ParameterSet&, FileName const&, int) const;
    int WriteParams(Cpptraj::Parm::ParameterSet&, FileName const&, int) const;
  private:
    int ReadInput(std::string&, BufferedLine&) const;
};
#endif
