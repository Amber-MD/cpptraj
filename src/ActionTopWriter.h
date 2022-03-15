#ifndef INC_ACTIONTOPWRITER_H
#define INC_ACTIONTOPWRITER_H
#include <string>
// Forward declares
class ActionSetup;
class ArgList;
class CoordinateInfo;
class Topology;
/// Class to hold common functionality for actions that will write modified topologies.
class ActionTopWriter {
  public:
    ActionTopWriter();
    ~ActionTopWriter();

    static const char* Keywords();
    static const char* Options();

    /// Parse arguments.
    int InitTopWriter(ArgList&, const char*, int);
    /// Write options to stdout.
    void PrintOptions() const;
    /// Write the Topology to file(s)
    int WriteTops(Topology const&) const;
    /// Remove box information from Topology/CoordinateInfo
    int ModifyActionState(ActionSetup&, Topology*);
  private:
    std::string prefix_;      ///< Prefix for writing topology as <prefix>.<originalname>
    std::string parmoutName_; ///< Output topology file name
    std::string parmOpts_;    ///< Topology file write args
    std::string typeStr_;     ///< Label for kind of topology being written.
    int debug_;               ///< Debug level to pass to Topology file writer.
    bool removeBoxInfo_;       ///< If true, remove box information from Topology/CoordinateInfo
    CoordinateInfo* newCinfo_; ///< Hold any modified CoordinateInfo

    static const char* keywords_;
    static const char* options_;
};
#endif
