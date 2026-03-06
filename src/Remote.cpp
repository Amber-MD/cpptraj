#include "Remote.h"
#include "CpptrajStdio.h"
#include "FileName.h"
#include <cstdlib> // system

using namespace Cpptraj;

/** The remote download command. */
std::string Remote::cmd_ = "";

/** Flag for setting output file for remote command. */
std::string Remote::oflag_ = "";

/** CONSTRUCTOR */
Remote::Remote() :
  overwrite_(false),
  debug_(0)
{
  setRemoteDownloadCommand();
}

/** CONSTRUCTOR - base url */
Remote::Remote(std::string const& baseUrl) :
  overwrite_(false),
  debug_(0),
  url_(baseUrl)
{
  setRemoteDownloadCommand();
}

/** Set whether to overwrite existing files or not. */
void Remote::SetOverwrite(bool b) {
  overwrite_ = b;
}

/** Set debug level */
void Remote::SetDebug(int d) {
  debug_ = d;
}

/** Set the remote download command. */
int Remote::setRemoteDownloadCommand() {
  if (!cmd_.empty()) return 0;
  // First try curl
# ifdef _WIN32
  int err = system("curl --version > NUL 2>&1");
# else
  int err = system("curl --version > /dev/null 2>&1");
# endif
  if (err == 0) {
    mprintf("\tcurl found.\n");
    cmd_.assign("curl -s --show-error -f -L ");
    oflag_.assign("-o ");
  } else {
#   ifdef _WIN32
    err = system("wget --version > NUL 2>&1");
#   else
    err = system("wget --version > /dev/null 2>&1");
#   endif
    if (err == 0) {
      mprintf("\twget found.\n");
      cmd_.assign("wget --quiet ");
      oflag_.assign("-O ");
    } else {
      mprinterr("Error: No working remote command found.\n");
      return 1;
    }
  }
  return 0;
}

/** Download file from base URL */
int Remote::DownloadFile(std::string const& fname, std::string const& outputFnameIn)
const
{
  if (cmd_.empty()) {
    mprinterr("Error: No remote download command set.\n");
    return 1;
  }
  if (fname.empty()) {
    mprinterr("Error: No file name to download specified.\n");
    return 1;
  }
  FileName fileName(fname);
  // Set output file name
  FileName outputFname;
  if (outputFnameIn.empty())
    outputFname = fileName.Base();
  else
    outputFname.SetFileName( outputFnameIn );
  bool fileExists = File::Exists( outputFname );
  if (overwrite_) {
    if (fileExists)
      mprintf("Warning: Overwriting existing file '%s'\n", outputFname.full());
  } else {
    if (fileExists) {
      mprintf("Warning: Not overwriting existing file '%s'\n", outputFname.full());
      return 0;
    }
  }
  // Set remote URL
  std::string remoteUrl = url_ + "/" + fileName.Full();
  mprintf("\t %s => %s\n", remoteUrl.c_str(), outputFname.full());
  // Download File
  std::string remoteCmd = cmd_ + oflag_ + outputFname.Full() + " " + remoteUrl;
  if (debug_ > 0)
    mprintf("DEBUG: %s\n", remoteCmd.c_str());
  int err = system(remoteCmd.c_str());
  if (err != 0) {
    mprinterr("Error: Could not download %s => %s\n", remoteUrl.c_str(), outputFname.full());
    // FIXME: wget will leave behind empty files here
    if (File::Exists(outputFname))
      File::Remove(outputFname);
    return 1;
  }

  return 0;
}

/** Download file from base URL */
int Remote::DownloadFile(std::string const& fname) const {
  return DownloadFile(fname, "");
}
