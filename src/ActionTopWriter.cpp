#include "ActionTopWriter.h"
#include "ArgList.h"
#include "CpptrajStdio.h"
#include "Topology.h"
#include "ParmFile.h"

static const char* keywords_ = "\t[outprefix <prefix>] [parmout <filename>]\n"
                               "\t[parmopts <comma-separated-list>]\n";
static const char* options_ = 
          "    outprefix <prefix> : Write modified topology to <prefix>.<originalname>\n"
          "    parmout <filename> : Write modified topology to <filename>\n"
          "    parmopts <list>    : Options for writing topology file\n";

const char* ActionTopWriter::Keywords() { return keywords_; }
const char* ActionTopWriter::Options() { return options_; }

/** CONSTRUCTOR */
ActionTopWriter::ActionTopWriter() :
  debug_(0)
{}

/** Parse arguments. */
int ActionTopWriter::InitTopWriter(ArgList& argIn, const char* typeStrIn, int debugIn) {
  debug_ = debugIn;
  prefix_ = argIn.GetStringKey("outprefix");
  parmoutName_ = argIn.GetStringKey("parmout");
  parmOpts_ = argIn.GetStringKey("parmopts");
  if (typeStrIn==0) {
    mprinterr("Internal Error: No type string given for ActionTopWriter::InitTopWriter\n");
    return 1;
  }
  typeStr_.assign(typeStrIn);
  return 0;
}

/** Print arguments to stdout. */
void ActionTopWriter::PrintOptions() const {
  if (!prefix_.empty())
    mprintf("\tWriting '%s' topology file with prefix '%s'\n", typeStr_.c_str(), prefix_.c_str());
  if (!parmoutName_.empty())
    mprintf("\tWriting '%s' topology file with name '%s'\n", typeStr_.c_str(), parmoutName_.c_str());
  if (!parmOpts_.empty())
    mprintf("\tUsing the following write options: %s\n", parmOpts_.c_str());
}

/** Write Topology to topology file(s).
  * \param topIn Topology to write.
  * \return 0 if OK, 1 if prefix write failed, 2 if parmout write failed, 3 if both failed.
  */
int ActionTopWriter::WriteTops(Topology const& topIn) const {
  ArgList writeOpts;
  int err = 0;
  if (!parmOpts_.empty())
    writeOpts.SetList(parmOpts_, ",");
  if (!prefix_.empty()) {
    ParmFile pfile;
    if ( pfile.WritePrefixTopology(topIn, prefix_, writeOpts, ParmFile::UNKNOWN_PARM, debug_) )
    {
      mprinterr("Error: Could not write out '%s' topology file with prefix '%s'.\n",
                typeStr_.c_str(), prefix_.c_str());
      err = 1;
    }
  }
  if (!parmoutName_.empty()) {
    ParmFile pfile;
    if (pfile.WriteTopology(topIn, parmoutName_, writeOpts, ParmFile::UNKNOWN_PARM, debug_))
    {
      mprinterr("Error: Could not write out '%s' topology file '%s'\n",
                typeStr_.c_str(), parmoutName_.c_str());
      err += 2;
    }
  }
  return err;
}
