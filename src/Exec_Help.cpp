#include "Exec_Help.h"
#include "CpptrajStdio.h"
#include "Command.h"

Exec_Help::Exec_Help() : Exec(GENERAL) { }

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

Exec::RetType Exec_Help::Execute(CpptrajState& State, ArgList& argIn) {
  ArgList arg = argIn;
  arg.RemoveFirstArg();
  if (arg.empty())
    // NONE in this context means list all commands
    Command::ListCommands( NONE );
  else {
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
    else
      cmd.Help();
  }
  return CpptrajState::OK;
}
