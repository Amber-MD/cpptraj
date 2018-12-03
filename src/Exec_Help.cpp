#include "Exec_Help.h"
#include "CpptrajStdio.h"
#include "Command.h"

void Exec_Help::Help() const {
  mprintf("\t[ { all |\n"
          "\t    <cmd> |\n"
          "\t    <command category> |\n"
          "\t    Formats [{trajin|trajout|readdata|writedata} [<format key>]]\n"
          "\t   } ]\n"
          "\tCommand Categories:");
  for (int i = 0; i != (int)DEPRECATED; i++) {
    const char* catKey = ObjKeyword((Otype)i);
    if (catKey != 0)
      mprintf(" %s", catKey);
  }
  mprintf("\n");
  mprintf("  all                : Print all known commands.\n"
          "  <cmd>              : Print help for command <cmd>.\n"
          "  <command category> : Print all commands in specified category.\n"
          "  Formats            : Help for file formats.\n");
}

/** Print help for file formats. */
int Exec_Help::Formats(ArgList& argIn) const {
  std::string ftype = argIn.GetStringNext();
  if (ftype == "trajin") {
    std::string fkey = argIn.GetStringNext();
    if (fkey.empty()) mprintf("    Available input trajectory formats:\n");
    TrajectoryFile::ReadOptions(fkey);
    if (fkey.empty())
      mprintf("    Use 'help Formats trajin <format key> for format-specific help.\n");
  } else if (ftype == "trajout") {
    std::string fkey = argIn.GetStringNext();
    if (fkey.empty()) mprintf("    Available output trajectory formats:\n");
    TrajectoryFile::WriteOptions(fkey);
    if (fkey.empty())
      mprintf("    Use 'help Formats trajout <format key> for format-specific help.\n");
  } else if (ftype == "readdata") {
    std::string fkey = argIn.GetStringNext();
    if (fkey.empty()) mprintf("    Available input datafile formats:\n");
    DataFile::ReadOptions(fkey);
    if (fkey.empty())
      mprintf("    Use 'help Formats readdata <format key> for format-specific help.\n");
  } else if (ftype == "writedata") {
    std::string fkey = argIn.GetStringNext();
    if (fkey.empty()) mprintf("    Available output datafile formats:\n");
    DataFile::WriteOptions(fkey);
    if (fkey.empty())
      mprintf("    Use 'help Formats writedata <format key> for format-specific help.\n");
  } else {
    mprintf("\t[{trajin|trajout|readdata|writedata} [<format keyword>]]\n"
            "  trajin    : List available input trajectory formats.\n"
            "  trajout   : List available output trajectory formats.\n"
            "  readdata  : List of available input datafile formats.\n"
            "  writedata : List of available output datafile formats.\n"
            "  <format keyword> : If specified provide specific help for that format.\n");
  }
  return 1;
}

/** \return 1 if a help topic was found, 0 otherwise. */
int Exec_Help::Topics(ArgList& argIn) const {
  // By convention, Topics will start with uppercase letters and
  // commands will start with lower case.
  if ( isupper(argIn[0][0]) ) {
    if (argIn[0].compare(0,6,"Format")==0)
      return Formats(argIn);
  }
  return 0;
}

Exec::RetType Exec_Help::Execute(CpptrajState& State, ArgList& argIn) {
  ArgList arg = argIn;
  arg.RemoveFirstArg();
  if (arg.empty())
    Help();
  else if (arg.CommandIs("all"))
    Command::ListCommands( NONE );
  else {
    arg.MarkArg(0);
    // Check for help topic.
    if (Topics(arg)) return CpptrajState::OK;
    // Check for command category.
    for (int i = 1; i != (int)DEPRECATED; i++) {
      if (arg.CommandIs( ObjKeyword((Otype)i) )) {
        Command::ListCommands( (Otype)i );
        return CpptrajState::OK;
      }
    }
    // Find help for specified command.
    Cmd const& cmd = Command::SearchToken( arg );
    if (cmd.Empty())
      mprintf("No help found for '%s'\n", arg.Command());
    else {
      if (cmd.Obj().Type() == DispatchObject::DEPRECATED)
        mprintf("Warning: '%s' is deprecated.\n", arg.Command());
      //arg.MarkArg(0);
      cmd.Help();
    }
  }
  return CpptrajState::OK;
}
