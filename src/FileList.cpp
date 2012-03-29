#include "FileList.h"

// CONSTRUCTOR
FileList::FileList() :
  debug_(0)
{}

// DESTRUCTOR
FileList::~FileList() {}

// FileList::AddNames()
void FileList::AddNames(char *filename, std::string &basename, std::string &tag) {
  fnames_.push_back( std::string(filename) );
  basenames_.push_back( basename );
  tags_.push_back( tag );
}

// FileList::FindName()
int FileList::FindName(char *nameIn) {
  std::string name(nameIn);
  return FindName(name);
}

// FileList::FindName()
int FileList::FindName(std::string &nameIn) {
  if (nameIn.empty()) return -1;
  // If first char of name is a bracket, assume tag.
  int idx = 0;
  if ( !tags_.empty() && nameIn[0]=='[' ) {
    for (std::vector<std::string>::iterator tag = tags_.begin();
                                          tag != tags_.end(); tag++)
    {
      if ( *tag == nameIn ) return idx;
      ++idx;
    }
  }
  // If not a tag, prefer full filename match over base filename
  idx = 0;
  if (!fnames_.empty()) {
    for (std::vector<std::string>::iterator name = fnames_.begin();
                                            name != fnames_.end(); name++)
    {
      if ( *name == nameIn ) return idx;
      ++idx;
    }
  }
  // Last, try Base filenames
  idx = 0;
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

