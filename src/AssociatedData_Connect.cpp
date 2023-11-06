#include "AssociatedData_Connect.h"
#include "CpptrajStdio.h"
#include "ArgList.h"

/** CONSTRUCTOR - head and tail atom */
AssociatedData_Connect::AssociatedData_Connect(int head, int tail) :
  AssociatedData(CONNECT)
{
  connect_.push_back( head );
  connect_.push_back( tail );
}

const char* AssociatedData_Connect::HelpText = "[head <head atom>] [tail <tail atom>]";

int AssociatedData_Connect::ProcessAdataArgs(ArgList& argIn) {
  // Expect user atom args to start from 1
  int head = argIn.getKeyInt("head", 0) - 1;
  int tail = argIn.getKeyInt("tail", 0) - 1;
  if (head < 0) head = -1;
  if (tail < 0) tail = -1;
  if (head == -1 && tail == -1) {
    mprinterr("Error: Either at least 'head' or 'tail' must be specified.\n");
    return 1;
  }
  connect_.push_back( head );
  connect_.push_back( tail );
  return 0;
}

void AssociatedData_Connect::Ainfo() const {
  if (!connect_.empty()) {
    mprintf(" (Connect Atoms:");
    for (Iarray::const_iterator it = connect_.begin(); it != connect_.end(); ++it)
      mprintf(" %i", *it + 1);
    mprintf(")");
  }
  return;
}
