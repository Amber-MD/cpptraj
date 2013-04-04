#include "FileList.h"

// CONSTRUCTOR
FileList::FileList() :
  debug_(0)
{}

// DESTRUCTOR
FileList::~FileList() {}

void FileList::Clear() {
  fnames_.clear();
  tags_.clear();
}

// FileList::AddNameWithTag()
void FileList::AddNameWithTag(FileName const& filename, std::string const& tag) 
{
  fnames_.push_back( filename );
  tags_.push_back( tag );
}

// FileList::AddFilename()
void FileList::AddFilename(FileName const& filename) {
  fnames_.push_back( filename );
}

// FileList::FindName()
int FileList::FindName(std::string const& nameIn) const {
  if (nameIn.empty()) return -1;
  // If first char of name is a bracket, assume tag.
  int idx = 0;
  if ( !tags_.empty() && nameIn[0]=='[' ) {
    for (std::vector<std::string>::const_iterator tag = tags_.begin();
                                                  tag != tags_.end(); ++tag)
    {
      if ( *tag == nameIn ) return idx;
      ++idx;
    }
  }
  if (!fnames_.empty()) {
    // If not a tag, prefer full filename match over base filename
    idx = 0;
    for (std::vector<FileName>::const_iterator name = fnames_.begin();
                                               name != fnames_.end(); ++name)
    {
      if ( (*name).Full() == nameIn ) return idx;
      ++idx;
    }
    // Last, try Base filenames. TODO: Move this check to previous loop
    idx = 0;
    for (std::vector<FileName>::const_iterator name = fnames_.begin();
                                               name != fnames_.end(); ++name)
    {
      if ( (*name).Base() == nameIn ) return idx;
      ++idx;
    }
  }
  return -1;
}
