#ifndef __WIN32
#   include <wordexp.h>
#endif
#include <cstdio>  // FILE, fopen
#include <cerrno>  // fileErrMsg, errno
#include <cstring> // fileErrMsg, strerror
#include "FileName.h"
#include "CpptrajStdio.h"

#ifndef __WIN32
static void WexpErr(int err) {
  switch ( err ) {
    case WRDE_BADCHAR :
      mprinterr("Error: Illegal occurrence of newline or one of |, &, ;, <, >, (, ), {, }.\n");
      break;
    //case WRDE_BADVAL
    case WRDE_CMDSUB :
      mprinterr("Error: Command substitution is not allowed in file names.\n");
      break;
    case WRDE_NOSPACE :
      mprinterr("Error: Out of memory.\n");
      break;
    case WRDE_SYNTAX :
      mprinterr("Error: Bad syntax (unbalanced parentheses, unmatched quotes.\n");
      break;
  }
}
#endif /* __WIN32 */

// COPY CONSTRUCTOR
FileName::FileName( const FileName& rhs ) : fullPathName_(rhs.fullPathName_),
  baseName_(rhs.baseName_), extension_(rhs.extension_),
  compressExt_(rhs.compressExt_), dirPrefix_(rhs.dirPrefix_) {}

// ASSIGNMENT
FileName& FileName::operator=(const FileName& rhs) {
  if (this != &rhs) {
    fullPathName_ = rhs.fullPathName_;
    baseName_ = rhs.baseName_;
    extension_ = rhs.extension_;
    compressExt_ = rhs.compressExt_;
    dirPrefix_ = rhs.dirPrefix_;
  }
  return *this;
}

// FileName::clear()
void FileName::clear() {
  fullPathName_.clear();
  baseName_.clear();
  extension_.clear();
  compressExt_.clear();
  dirPrefix_.clear();
}

bool FileName::MatchFullOrBase(std::string const& rhs) const {
  if (!fullPathName_.empty()) {
    // Prefer full filename match.
    if (fullPathName_ == rhs) return true;
    if (baseName_     == rhs) return true;
  }
  return false;
}

int FileName::SetFileName(std::string const& nameIn) {
  // null filename allowed (indicates STDIN/STDOUT)
  if (nameIn.empty()) {
    clear();
    return 0;
  }
#ifndef __WIN32
  wordexp_t expanded;
  int err = wordexp( nameIn.c_str(), &expanded, WRDE_NOCMD );
  WexpErr( err );
  if (err == 0) {
    if (expanded.we_wordc < 1) { // Sanity check
      mprinterr("Internal Error: Word expansion failed.\n");
      err = 1;
    } else
      err = SetFileName_NoExpansion( expanded.we_wordv[0] );
    wordfree( &expanded );
  }
  return err;
#else
  SetFileName_NoExpansion(nameIn);
  return 0;
#endif
}

int FileName::SetFileName_NoExpansion(std::string const& nameIn) {
  // null filename allowed (indicates STDIN/STDOUT)
  if (nameIn.empty()) {
    clear();
    return 0;
  }
  // Assign filename with full path
  fullPathName_.assign( nameIn );
  // Get position of last occurence of '/' to determine base filename
  size_t found = fullPathName_.find_last_of("/");
  if (found == std::string::npos) {
    baseName_ = fullPathName_;
    dirPrefix_.clear();
  } else {
    baseName_ = fullPathName_.substr(found+1);
    dirPrefix_ = fullPathName_.substr(0, found+1);
  }
  // Get the filename extension
  found = baseName_.find_last_of(".");
  if (found == std::string::npos) {
    extension_.clear();
  } else {
    extension_ = baseName_.substr(found);
  }
  // See if the extension is one of the 2 recognized compression extensions.
  // If file has a compression format extension, look for another extension.
  if ( extension_ == ".gz" || extension_ == ".bz2" ) {
    compressExt_ = extension_;
    // Get everything before the compressed extension
    std::string compressPrefix = baseName_.substr(0,found);
    // See if there is another extension
    found = compressPrefix.find_last_of(".");
    if (found == std::string::npos)
      // No other extension
      extension_.clear();
    else
      extension_ = compressPrefix.substr(found);
  } else
    compressExt_.clear();
  return 0;
}

int FileName::Append( std::string const& suffix ) {
  if (fullPathName_.empty()) return 1;
  fullPathName_.append( suffix );
  baseName_.append( suffix );
  return 0;
}

FileName FileName::AppendFileName( std::string const& suffix ) const {
  FileName out( *this );
  out.Append( suffix );
  return out;
}

//TODO make this more efficient by just modifying full and base names
FileName FileName::PrependFileName( std::string const& prefix ) const {
  FileName out;
  out.SetFileName_NoExpansion(dirPrefix_ + prefix + baseName_);
  return out;
}

FileName FileName::PrependExt( std::string const& extPrefix ) const {
  FileName out( *this );
  // Find location of extension.
  size_t found = out.baseName_.rfind( extension_ );
  // Remove extension.
  out.baseName_.resize( found );
  // Insert extPrefix to just before extension and re-add extension.
  out.baseName_.append( extPrefix + extension_ + compressExt_ );
  // Update full path name.
  out.fullPathName_ = dirPrefix_ + out.baseName_;
  //mprintf("DEBUG: fullPathName= '%s'\n"
  //        "       baseName=     '%s'\n"
  //        "       extension=    '%s'\n"
  //        "       compressExt=  '%s'\n"
  //        "       dirPrefix=    '%s'\n",
  //        out.fullPathName_.c_str(), out.baseName_.c_str(), out.extension_.c_str(),
  //        out.compressExt_.c_str(), out.dirPrefix_.c_str());
  return out;
}

// =============================================================================
File::NameArray File::ExpandToFilenames(std::string const& fnameArg) {
  NameArray fnames;
#ifdef __WIN32
  fnames.push_back( fnameArg );
#else
  if (fnameArg.empty()) return fnames;
  wordexp_t expanded;
  int err = wordexp( fnameArg.c_str(), &expanded, WRDE_NOCMD );
  WexpErr( err );
  if ( err == 0 ) {
    for (unsigned int i = 0; i != expanded.we_wordc; i++) {
      if (expanded.we_wordv[i] == 0)
        mprinterr("Internal Error: Bad expansion at %i\n", i);
      else {
        FileName fn;
        fn.SetFileName_NoExpansion( expanded.we_wordv[i] );
        fnames.push_back( fn );
      }
    }
    wordfree( &expanded );
  }
#endif /* __WIN32 */
  return fnames;
}

static std::string fileErrMsg_ = std::string("");

void File::ErrorMsg(const char* fname) {
  mprinterr("Error: '%s': %s\n", fname, fileErrMsg_.c_str());
}

bool File::Exists(FileName const& fn) {
  if (!fn.empty()) {
    FILE* infile = fopen(fn.full(), "rb");
    if (infile==0) {
      fileErrMsg_.assign( strerror(errno) );
      return false;
    }
    fclose(infile);
    return true;
  }
  return false;
}

bool File::Exists(std::string const& fname) {
  return File::Exists( FileName(fname) );
}
