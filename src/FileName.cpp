#include <wordexp.h>
#include <cstdio>
#include <cstdlib>
#include <cstring> // strerror
#include <cerrno>  // errno
//#ifndef __PGI
//#  include <glob.h>  // For tilde expansion
//#endif
#include "FileName.h"
#include "CpptrajStdio.h"

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
/*# ifdef __PGI
  // NOTE: It seems some PGI compilers do not function correctly when glob.h
  //       is included and large file flags are set. Try to expand any tildes.
  char* HOME = getenv("HOME");
  std::string fname;
  if (HOME != 0) {
    for (std::string::const_iterator c = nameIn.begin(); c != nameIn.end(); ++c) {
      if (*c == '~')
        fname.append( HOME );
      else
        fname += *c;
    }
  } else
    fname = nameIn;
  return SetFileName_NoExpansion( fname );
//# else*/
/*
  globbuf.gl_offs = 1;
  std::string fname;
  int err = glob(nameIn.c_str(), GLOB_TILDE, NULL, &globbuf);
  if ( err == GLOB_NOMATCH )
    //mprinterr("Error: '%s' does not exist.\n", filenameIn); // Make silent
    return returnFilename;
  else if ( err != 0 )
    mprinterr("Internal Error: glob() failed for '%s' (%i)\n", filenameIn.c_str(), err);
  else {
    returnFilename.assign( globbuf.gl_pathv[0] );
    globfree(&globbuf);
  }
  return returnFilename;
//# endif
*/
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

int FileName::AppendFileName( std::string const& suffix ) {
  if (fullPathName_.empty()) return 1;
  fullPathName_.append( suffix );
  baseName_.append( suffix );
}

// =============================================================================
File::NameArray File::ExpandToFilenames(std::string const& fnameArg) {
  NameArray fnames;
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
  return fnames;
/*
  StrArray fnames;
  if (fnameArg.empty()) return fnames;
# ifdef __PGI
  // NOTE: It seems some PGI compilers do not function correctly when glob.h
  //       is included and large file flags are set. Just disable globbing
  //       for PGI and return a copy of filenameIn.
  // Check for any wildcards in fnameArg
  if ( fnameArg.find_first_of("*?[]") != std::string::npos )
    fprintf(stdout,"Warning: Currently wildcards in filenames not supported with PGI compilers.\n");
  fnames.push_back( fnameArg );
# else
  glob_t globbuf;
  int err = glob(fnameArg.c_str(), GLOB_TILDE, NULL, &globbuf );
  //printf("DEBUG: %s matches %zu files.\n", fnameArg.c_str(), (size_t)globbuf.gl_pathc);
  if ( err == 0 ) {
    for (unsigned int i = 0; i < (size_t)globbuf.gl_pathc; i++)
      fnames.push_back( globbuf.gl_pathv[i] );
  } else if (err == GLOB_NOMATCH )
    fprintf(stderr,"Error: %s matches no files.\n", fnameArg.c_str());
  else
    fprintf(stderr,"Error: occurred trying to find %s\n", fnameArg.c_str());
  if ( globbuf.gl_pathc > 0 ) globfree(&globbuf);
# endif
  return fnames;
*/
}

bool File::Exists(std::string const& fname) {
  FileName fn( fname );
  if (!fn.empty()) {
    FILE* infile = fopen(fn.full(), "rb");
    if (infile==0) {
      mprinterr("Error: File '%s': %s\n", fname.c_str(), strerror( errno ));
      return false;
    }
    fclose(infile);
    return true;
  }
  return false;
}
