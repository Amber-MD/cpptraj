#ifndef INC_COMMAND_H
#define INC_COMMAND_H
#include "CmdList.h"
#include "CpptrajState.h"
class ControlBlock;
class Command {
  public:
    static void Init();
    static void Free();
    /// List commands of given type, or all if type is NONE
    static void ListCommands(DispatchObject::Otype);
    /// \return command corresponding to first argument in ArgList.
    static Cmd const& SearchToken(ArgList&);
    /// \return command of given type corresponding to given command key. 
    static Cmd const& SearchTokenType(DispatchObject::Otype, const char*);
    /// \return true if unterminated control block(s) exist.
    static bool UnterminatedControl();
    /// Execute command, modifies given CpptrajState
    static CpptrajState::RetType Dispatch(CpptrajState&, std::string const&);
    /// Read input commands from given file, modifies given CpptrajState.
    static CpptrajState::RetType ProcessInput(CpptrajState&, std::string const&);
    /// \return Pointer to command name address.
    static const char* CmdToken(int idx) { return names_[idx]; }
  private:
    /// Add a command to the list of commands
    static void AddCmd(DispatchObject*, Cmd::DestType, int, ...);
    /// Search for command of given type with given name; can be silent.
    static Cmd const& SearchTokenType(DispatchObject::Otype, const char*, bool);
    /// List all commands for given type
    static void ListCommandsForType(DispatchObject::Otype);
    /// Clear all existing control blocks
    static void ClearControlBlocks();
    /// Add a new control block
    static int AddControlBlock(ControlBlock*, CpptrajState&, ArgList&);
    /// Execute specified control block
    static int ExecuteControlBlock(int, CpptrajState&);
    /// Execute given command
    static CpptrajState::RetType ExecuteCommand(CpptrajState&, ArgList const&);

    static CmdList commands_; ///< Master list of commands.
    static const Cmd EMPTY_;  ///< Empty command.
    typedef std::vector<const char*> Carray;
    static Carray names_;     ///< Array of pointers to all command names, for ReadLine
    typedef std::vector<ControlBlock*> CtlArray;
    static CtlArray control_; ///< Array of control blocks
    static int ctlidx_;       ///< Point to current control block
};
#endif
