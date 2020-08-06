#ifndef INC_ACTIONTOPWRITER_H
#define INC_ACTIONTOPWRITER_H
#include <string>
// Forward declares
class ArgList;
/// Class to hold common functionality for actions that will write modified topologies.
class ActionTopWriter {
  public:
    ActionTopWriter();

    static const char* Keywords();
    static const char* Options();

    /// Parse arguments.
    int InitTopWriter(ArgList&);
    /// Write options to stdout.
    void PrintOptions(const char*) const;
  private:
    std::string prefix_;      ///< Prefix for writing topology as <prefix>.<originalname>
    std::string parmoutName_; ///< Output topology file name
    std::string parmOpts_;    ///< Topology file write args
};
#endif
