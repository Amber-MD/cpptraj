#include "Exec_Commands.h"
#include "CpptrajStdio.h"

void Exec_Run::Help() const {
  mprintf("  Process all trajectories currently in input trajectory list.\n"
          "  All actions in action list will be run on each frame.\n"
          "  If not processing ensemble input, all analyses in analysis\n"
          "  list will be run after trajectory processing.\n");
}
