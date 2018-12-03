#include "Exec_Help.h"
#include "CpptrajStdio.h"
#include "Command.h"

void Exec_Help::Help() const {
  mprintf("\t[{ <cmd> | <category>}]\n\tCategories:");
  for (int i = 0; i != (int)DEPRECATED; i++) {
    const char* catKey = ObjKeyword((Otype)i);
    if (catKey != 0)
      mprintf(" %s", catKey);
  }
  mprintf("\n");
  mprintf("  With no arguments list all known commands, otherwise display help for specified\n"
          "  command. If a category is specified list only commands in that category.\n");
}

/** Print help for file formats. */
int Exec_Help::Formats(ArgList& argIn) const {
  std::string ftype = argIn.GetStringNext();
  if (ftype == "trajin") {
    mprintf("    Available input trajectory formats:\n");
    TrajectoryFile::ReadOptions(argIn.GetStringNext());
  } else if (ftype == "trajout") {
    mprintf("    Available output trajectory formats:\n");
    TrajectoryFile::WriteOptions(argIn.GetStringNext());
  } else {
    mprintf("\t[{trajin|trajout} [<format keyword>]]\n"
            "  trajin           : List available input trajectory formats.\n"
            "  trajout          : List available output trajectory formats.\n"
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
    // NONE in this context means list all commands
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
      cmd.Help( arg );
    }
  }
  return CpptrajState::OK;
}
