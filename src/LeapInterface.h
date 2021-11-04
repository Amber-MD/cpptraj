#ifndef INC_LEAPINTERFACE_H
#define INC_LEAPINTERFACE_H
#include <vector>
#include <string>
class FileName;
namespace Cpptraj {
/// Provide an interface to generate a topology from Amber LEaP
class LeapInterface {
  public:
    /// CONSTRUCTOR
    LeapInterface();
    /// Add input file
    int AddInputFile(std::string const&);
    /// Clear input files
    void ClearInputFiles();
    /// Run leap
    int RunLeap() const;
  private:
    typedef std::vector<std::string> Sarray;

    /// Execute leap process
    int execute_leap(FileName const&) const;

    Sarray input_files_; ///< Array of leap input files to source
    std::string leapOutName_; ///< Leap output file name
};
}
#endif
