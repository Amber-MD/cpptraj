#include "Exec_Help.h"
#include "CpptrajStdio.h"
#include "Command.h"
#include "ParmFile.h"

void Exec_Help::Help() const {
  mprintf("\t[ { all |\n"
          "\t    <cmd> |\n"
          "\t    <command category> |\n"
          "\t    Formats [{read|write}] |\n"
          "\t    Formats [{trajin|trajout|readdata|writedata|parm|parmwrite} [<fmt key>]] |\n"
          "\t    Masks\n"
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
          "  Formats            : Help for file formats.\n"
          "  Masks              : Help for mask syntax.\n");
}

/** Print help for file formats. */
int Exec_Help::Formats(ArgList& argIn) const {
  std::string ftype = argIn.GetStringNext();
  if (ftype == "trajin") {
    std::string fkey = argIn.GetStringNext();
    if (fkey.empty()) mprintf("    *** Available input trajectory formats ***\n");
    TrajectoryFile::ReadOptions(fkey);
    if (fkey.empty())
      mprintf("    Use 'help Formats trajin <format key> for format-specific help.\n");
  } else if (ftype == "trajout") {
    std::string fkey = argIn.GetStringNext();
    if (fkey.empty()) mprintf("    *** Available output trajectory formats ***\n");
    TrajectoryFile::WriteOptions(fkey);
    if (fkey.empty())
      mprintf("    Use 'help Formats trajout <format key> for format-specific help.\n");
  } else if (ftype == "readdata") {
    std::string fkey = argIn.GetStringNext();
    if (fkey.empty()) mprintf("    *** Available input datafile formats ***\n");
    DataFile::ReadOptions(fkey);
    if (fkey.empty())
      mprintf("    Use 'help Formats readdata <format key> for format-specific help.\n");
  } else if (ftype == "writedata") {
    std::string fkey = argIn.GetStringNext();
    if (fkey.empty()) mprintf("    *** Available output datafile formats ***\n");
    DataFile::WriteOptions(fkey);
    if (fkey.empty())
      mprintf("    Use 'help Formats writedata <format key> for format-specific help.\n");
  } else if (ftype == "parm") {
    std::string fkey = argIn.GetStringNext();
    if (fkey.empty()) mprintf("    *** Available input topology formats ***\n");
    ParmFile::ReadOptions(fkey);
    if (fkey.empty())
      mprintf("    Use 'help Formats parm <format key> for format-specific help.\n");
  } else if (ftype == "parmwrite") {
    std::string fkey = argIn.GetStringNext();
    if (fkey.empty()) mprintf("    *** Available output topology formats ***\n");
    ParmFile::WriteOptions(fkey);
    if (fkey.empty())
      mprintf("    Use 'help Formats parmwrite <format key> for format-specific help.\n");
  } else if (ftype == "read") {
    mprintf("    *** Available input trajectory formats ***\n");
    TrajectoryFile::ReadOptions("");
    mprintf("    *** Available input datafile formats ***\n");
    DataFile::ReadOptions("");
    mprintf("    *** Available input topology formats ***\n");
    ParmFile::ReadOptions("");
  } else if (ftype == "write") {
    mprintf("    *** Available output trajectory formats ***\n");
    TrajectoryFile::WriteOptions("");
    mprintf("    *** Available output datafile formats ***\n");
    DataFile::WriteOptions("");
    mprintf("    *** Available output topology formats ***\n");
    ParmFile::WriteOptions("");
  } else {
    mprintf("\t[{read|write}] |\n"
            "\t[{trajin|trajout|readdata|writedata|parm|parmwrite} [<fmt key>]]\n"
            "  read      : List available input formats.\n"
            "  write     : List available output formats.\n"
            "  trajin    : List available input trajectory formats.\n"
            "  trajout   : List available output trajectory formats.\n"
            "  readdata  : List available input datafile formats.\n"
            "  writedata : List available output datafile formats.\n"
            "  parm      : List available input topology formats.\n"
            "  parmwrite : List available output topology formats.\n"
            "  <fmt key> : If specified provide specific help for that format.\n");
  }
  return 1;
}

/** Help for atom masks. */
int Exec_Help::Masks(ArgList& argIn) const {
  mprintf("    CPPTRAJ mask syntax.\n"
          "  *** Basic selection ***\n"
          "    @{list}  : Select atoms by number/name. E.g. '@1-5,12-17,20', '@CA', '@CA,C,O,N,H'\n"
          "    @%{list} : Select atom types. E.g. '@%%CT'\n"
          "    @/{list} : Select atom elements. E.g. '@/N'\n"
          "    :{list}  : Select residues by number/name. E.g. ':1-10,15,19-22', ':LYS', ':ASP,ALA'\n"
          "    :/{list} : Select residues by chain ID. E.g. ':/B', ':/A,D'.\n"
          "    :;{list} : Select by PDB residue number.\n"
          "    ^{list}  : Select by molecule number. E.g. '^1-10', '^2-4,8'\n"
          "  Combinations of atom/residue/molecule masks are interpreted as if 'AND'\n"
          "  is specified, e.g. ':WAT@O' is 'residues named WAT and atoms named O'.\n"
          "  *** Distance-based masks ***\n"
          "    <mask><distance op><distance>\n"
          "      <mask>        : Specify atoms to select around.\n"
          "      <distance op> : Distance operator.\n"
          "                      @< means 'atoms within'\n"
          "                      @> means 'atoms outside of'\n"
          "                      :< means 'residues within'\n"
          "                      :> means 'residues outside of'\n"
          "                      ^< means 'molecules within'\n"
          "                      ^> means 'molecules outside of'\n"
          "      <distance>    : Cutoff for distance operator.\n"
          "    E.g. ':11-17<@2.4' means 'select atoms within 2.4 Ang. distance of atoms\n"
          "      selected by ':11-17' (residues numbered 11 through 17).\n"
          "  *** Operators ***\n"
          "    ( ) : Open/close parentheses.\n"
          "    &   : AND operator.\n"
          "    |   : OR operator.\n"
          "    !   : NOT operator.\n"
          "  *** Wildcards ***\n"
          "    * : Zero or more characters. Can be used to select all. E.g. '@H*'.\n"
          "    = : Same as '*'\n"
          "    ? : Single character. E.g. ':?0', ':AS?'\n");
  return 1;
}

/** \return 1 if a help topic was found, 0 otherwise. */
int Exec_Help::Topics(ArgList& argIn) const {
  // By convention, Topics will start with uppercase letters and
  // commands will start with lower case.
  if ( isupper(argIn[0][0]) ) {
    if (argIn[0].compare(0,6,"Format")==0)
      return Formats(argIn);
    else if (argIn[0].compare(0,4,"Mask")==0)
      return Masks(argIn);
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
