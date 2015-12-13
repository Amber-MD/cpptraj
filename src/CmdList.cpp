#include "CmdList.h"

CmdList::~CmdList() { Clear(); }

void CmdList::Clear() {
  for (Carray::iterator it = cList_.begin(); it != cList_.end(); ++it)
    it->Clear();
  cList_.clear();
} 
