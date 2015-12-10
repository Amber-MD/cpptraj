#include "Cmd.h"

void Cmd::Clear() { if (object_ != 0) delete object_; }
