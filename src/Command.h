#ifndef INC_COMMAND_H
#define INC_COMMAND_H
#include "CmdList.h"
#include "CpptrajState.h"
#include "Control.h"
class Command {
  public:
    static void Init();
    static void Free();
    /// List commands of given type, or all if type is NONE
    static void ListCommands(DispatchObject::Otype);
    /// Add a command to the list of commands
    static void AddCmd(DispatchObject*, Cmd::DestType, int, ...);
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
    static Cmd const& SearchTokenType(DispatchObject::Otype, const char*, bool);
    static void ListCommandsForType(DispatchObject::Otype);
    static int AddControlBlock(Control*, CpptrajState&, ArgList&);
    static int ExecuteControlBlock(int, CpptrajState&, Control::Varray);
    static CpptrajState::RetType Dispatch(CpptrajState&, ArgList&);

    static CmdList commands_; ///< Master list of commands.
    static const Cmd EMPTY_;  ///< Empty command.
    typedef std::vector<const char*> Carray;
    static Carray names_; ///< Array of pointers to all command names, for ReadLine
    typedef std::vector<Control*> CtlArray;
    static CtlArray control_;
    static int ctlidx_;
};
#endif
