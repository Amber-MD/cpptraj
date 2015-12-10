#include "Exec_ReadInput.h"
#include "CpptrajStdio.h"
#include "Command.h"

void Exec_ReadInput::Help() const {
  mprintf("\t<filename>\n"
          "  Read commands from input file <filename>\n");
}

Exec::RetType Exec_ReadInput::Execute(CpptrajState& State, ArgList& argIn)
{
  // Next arg should be a filename. Not allowed to be blank in command.
  std::string inputFilename = argIn.GetStringNext();
  if (inputFilename.empty()) {
    mprinterr("Error: No input filename given.\n");
    return CpptrajState::ERR;
  }
  return Command::ProcessInput( State, inputFilename );
}
