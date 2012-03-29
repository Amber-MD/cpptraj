#include "FileList.h"

// CONSTRUCTOR
FileList::FileList() :
  debug_(0)
{}

// DESTRUCTOR
FileList::~FileList() {}

// FileList::FindTag()
int FileList::FindTag(std::string &tagIn) {
  if (tags_.empty()) return -1;
  int idx = 0;
  for (std::vector<std::string>::iterator tag = tags_.begin();
                                          tag != tags_.end(); tag++)
  {
    if ( *tag == tagIn ) return idx;
    ++idx;
  }
  return -1;
}

// FileList::FindName()
int FileList::FindName(std::string &nameIn) {
  // Prefer full filename match
  int idx = 0;
  if (!fnames_.empty()) {
    for (std::vector<std::string>::iterator name = fnames_.begin();
                                            name != fnames_.end(); name++)
    {
      if ( *name == nameIn ) return idx;
      ++idx;
    }
  }
  // Base filenames
  if (!basenames_.empty()) {
    for (std::vector<std::string>::iterator name = basenames_.begin();
                                            name != basenames_.end(); name++)
    {
      if ( *name == nameIn ) return idx;
      ++idx;
    }
  }
  return -1;
}

