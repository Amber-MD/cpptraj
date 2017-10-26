#ifndef _WIN32
#   include <wordexp.h>
#endif
#include <cstdio>  // FILE, fopen
#include <cerrno>  // fileErrMsg, errno
#include <cstring> // fileErrMsg, strerror
#include "FileName.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // validInteger, convertToInteger, integerToString

#ifndef _WIN32
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
#endif /* _WIN32 */

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
#ifndef _WIN32
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
// ----- RepName Class ---------------------------------------------------------
// RepName CONSTRUCTOR
File::RepName::RepName(FileName const& fname, int debugIn) :
  extChar_('.')
{
  if (debugIn > 1)
    mprintf("\tREMDTRAJ: FileName=[%s]\n", fname.full());
  if ( fname.Ext().empty() ) {
    mprinterr("Error: Traj %s has no numerical extension, required for automatic\n"
              "Error:   detection of replica trajectories. Expected filename format is\n"
              "Error:   <Prefix>.<#> (with optional compression extension), examples:\n"
              "Error:   Rep.traj.nc.000,  remd.x.01.gz etc.\n", fname.base());
    return;
  }
  // Split off everything before replica extension
  size_t found = fname.Full().rfind( fname.Ext() );
  Prefix_.assign( fname.Full().substr(0, found) );
  ReplicaExt_.assign( fname.Ext() ); // This should be the numeric extension
  // Remove leading '.'
  if (ReplicaExt_[0] == '.') ReplicaExt_.erase(0,1);
  CompressExt_.assign( fname.Compress() );
  if (debugIn > 1) {
    mprintf("\tREMDTRAJ: Prefix=[%s], #Ext=[%s], CompressExt=[%s]\n",
            Prefix_.c_str(), ReplicaExt_.c_str(), CompressExt_.c_str());
  }
  // CHARMM replica numbers are format <name>_<num>
  if ( !validInteger(ReplicaExt_) ) {
    size_t uscore = fname.Full().rfind('_');
    if (uscore != std::string::npos) {
      Prefix_.assign( fname.Full().substr(0, uscore) );
      ReplicaExt_.assign( fname.Full().substr(uscore+1) );
      extChar_ = '_';
      if (debugIn > 0)
        mprintf("\tREMDTRAJ: CHARMM style replica names detected, prefix='%s' ext='%s'\n",
                Prefix_.c_str(), ReplicaExt_.c_str());
    }
  }
  // Check that the numerical extension is valid.
  if ( !validInteger(ReplicaExt_) ) {
    mprinterr("Error: Replica extension [%s] is not an integer.\n", ReplicaExt_.c_str());
    Prefix_.clear(); // Empty Prefix_ indicates error.
    return;
  }
  ExtWidth_ = (int)ReplicaExt_.size();
  if (debugIn > 1)
    mprintf("\tREMDTRAJ: Numerical Extension width=%i\n", ExtWidth_);
  // Store lowest replica number
  lowestRepnum_ = convertToInteger( ReplicaExt_ );
  // TODO: Do not allow negative replica numbers?
  if (debugIn > 1)
    mprintf("\tREMDTRAJ: index of first replica = %i\n", lowestRepnum_);
}

/** \return Replica file name for given offset from lowest replica number. */
FileName File::RepName::RepFilename(int offset) const {
  FileName trajFilename;
  trajFilename.SetFileName_NoExpansion( Prefix_ + extChar_ +
                                        integerToString(lowestRepnum_ + offset, ExtWidth_) +
                                        CompressExt_ );
  return trajFilename;
}

// -----------------------------------------------------------------------------
File::NameArray File::ExpandToFilenames(std::string const& fnameArg) {
  NameArray fnames;
#ifdef _WIN32
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
#endif /* _WIN32 */
  return fnames;
}

/** Assuming lowest replica filename has been set, search for all other 
  * replica names assuming a naming scheme of '<PREFIX>.<EXT>[.<CEXT>]', 
  * where <EXT> is a numerical extension and <CEXT> is an optional 
  * compression extension. 
  * \return Found replica filenames, or an empty list on error. 
  */
File::NameArray File::SearchForReplicas(FileName const& fname, int debug) {
  NameArray replica_filenames;
  RepName repName(fname, debug);
  if (repName.Error()) return replica_filenames;
  // Search for a replica number lower than this. Correct functioning
  // of the replica code requires the file specified by trajin be the
  // lowest # replica.
  if (File::Exists( repName.RepFilename( -1 ) )) {
    mprintf("Warning: Replica# found lower than file specified with trajin.\n"
            "Warning:   Found \"%s\"; 'trajin remdtraj' requires lowest # replica.\n",
            repName.RepFilename( -1 ).full());
  }
  // Add lowest replica filename, search for and add all replicas higher than it.
  replica_filenames.push_back( fname );
  int rep_offset = 0;
  bool search_for_files = true;
  FileName trajFilename;
  while (search_for_files) {
    ++rep_offset;
    trajFilename = repName.RepFilename( rep_offset );
    //mprintf("\t\tChecking for %s\n", trajFilename.full());
    if (File::Exists( trajFilename ))
      replica_filenames.push_back( trajFilename );
    else
      search_for_files = false;
  }
  return replica_filenames;
}

/** Each rank searches for replica based on lowest replica number. */
File::NameArray File::SearchForReplicas(FileName const& fname, bool trajCommMaster,
                                        int ensRank, int ensSize, int debug)
{
  NameArray replica_filenames;
  RepName repName(fname, debug);
  if (repName.Error()) return replica_filenames;
  // TODO check for lower replica number?
  FileName replicaFilename = repName.RepFilename( ensRank );
  // Only traj comm masters actually check for files.
  if (trajCommMaster) {
    if (!File::Exists( replicaFilename )) {
      File::ErrorMsg( replicaFilename.full() );
      rprinterr("Error: File '%s' not accessible.\n", replicaFilename.full());
      return replica_filenames;
    }
  }
  // At this point each rank has found its replica. Populate filename array.
  for (int offset = 0; offset < ensSize; ++offset)
    replica_filenames.push_back( repName.RepFilename( offset ) );
  return replica_filenames;
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
