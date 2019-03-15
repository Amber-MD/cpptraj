#include "Action.h"
#include "CpptrajStdio.h"

void Action::CheckImageRotationWarning(ActionSetup const& setup, const char* desc) {
  if (setup.CoordInfo().TrajBox().Type() != Box::NOBOX) {
    mprintf("Warning: Coordinates are being rotated and box coordinates are present.\n"
            "Warning: Unit cell vectors are NOT rotated; imaging will not be possible\n"
            "Warning:  after %s is performed.\n");
  }
}
