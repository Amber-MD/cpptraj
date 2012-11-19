#include "FileList.h"

// CONSTRUCTOR
FileList::FileList() :
  debug_(0)
{}

// DESTRUCTOR
FileList::~FileList() {}

void FileList::Clear() {
  fnames_.clear();
  basenames_.clear();
  tags_.clear();
}

// FileList::AddNameWithTag()
void FileList::AddNameWithTag(std::string const& filename, 
                              std::string const& basename, 
                              std::string const& tag) 
{
  fnames_.push_back( filename );
  basenames_.push_back( basename );
  tags_.push_back( tag );
}

// FileList::AddFilename()
void FileList::AddFilename(std::string const& filename) {
  fnames_.push_back( filename );
}

// FileList::FindName()
int FileList::FindName(std::string const& nameIn) {
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

// FileList::HasNames()
bool FileList::HasNames() {
  if (!tags_.empty() || !fnames_.empty() || !basenames_.empty())
    return true;
  return false;
}
