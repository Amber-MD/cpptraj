#include "ActionTopWriter.h"
#include "ArgList.h"
#include "CpptrajStdio.h"

static const char* keywords_ = "\t[outprefix <prefix>] [parmout <filename>]\n"
                               "\t[parmopts <comma-separated-list>]\n";
static const char* options_ = 
          "    outprefix <prefix> : Write re-ordered topology to <prefix>.<originalname>\n"
          "    parmout <filename> : Write re-ordered topology to <filename>\n"
          "    parmopts <list>    : Options for writing topology file\n";

const char* ActionTopWriter::Keywords() { return keywords_; }
const char* ActionTopWriter::Options() { return options_; }

/** CONSTRUCTOR */
ActionTopWriter::ActionTopWriter() {}

/** Parse arguments. */
int ActionTopWriter::InitTopWriter(ArgList& argIn) {
  prefix_ = argIn.GetStringKey("outprefix");
  parmoutName_ = argIn.GetStringKey("parmout");
  parmOpts_ = argIn.GetStringKey("parmopts");
  return 0;
}

/** Print arguments to stdout. */
void ActionTopWriter::PrintOptions(const char* typeStr) const {
  if (!prefix_.empty())
    mprintf("\t%s will be output with prefix '%s'\n", typeStr, prefix_.c_str());
  if (!parmoutName_.empty())
    mprintf("\t%s will be output with name '%s'\n", typeStr, parmoutName_.c_str());
  if (!parmOpts_.empty())
    mprintf("\tUsing the following write options: %s\n", parmOpts_.c_str());
}
