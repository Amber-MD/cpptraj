#ifndef INC_COMMAND_H
#define INC_COMMAND_H
#include "Cmd.h"
/// This is a static class that determines how commands are handled.
/** To add a new Action/Analysis command, add the appropriate '#include'
  * to the top of Command.cpp and add an entry to Commands[] (search for 
  * INC_ACTION/INC_ANALYSIS respectively).
  */
class Command {
  public:
    static void ListCommands(Cmd::Ctype);
    static Cmd::TokenPtr SearchTokenType(Cmd::Ctype, ArgList const& argIn);
    static Cmd::TokenPtr SearchToken(ArgList&);
    static Cmd::RetType Dispatch(CpptrajState&, std::string const&);
    static Cmd::RetType ProcessInput(CpptrajState&, std::string const&);
    static Cmd::Token const& CmdToken(int idx)       { return Commands[idx]; }
    static const char* CommandCategoryKeyword(int i) { return CommandTitle[i]; }
  private:
    static void WarnDeprecated(Cmd::TokenPtr);
    static const char* CommandTitle[];
    /// Master list of commands.
    static const Cmd::Token Commands[];
};
#endif
