#include "LeapInterface.h"
#include "ArgList.h"
#include "CpptrajFile.h"
#include "CpptrajStdio.h"
#include "File_TempName.h"
#include "StringRoutines.h"
#include <cstdio>

using namespace Cpptraj;
/** CONSTRUCTOR */
LeapInterface::LeapInterface() : debug_(0) {}

/** CONSTRUCTOR - set debug level */
LeapInterface::LeapInterface(int d) : debug_(d) {}

/** Set up leap to run. */
int LeapInterface::AddInputFile(std::string const& fname) {
  if (fname.empty()) return 1;
  input_files_.push_back( fname );
  return 0;
}

/** Add command to run after input is sourced. */
int LeapInterface::AddCommand(std::string const& cmd) {
  if (cmd.empty()) return 1;
  commands_.push_back( cmd );
  return 0;
}

/** Clear all input files. */
void LeapInterface::ClearInputFiles() {
  input_files_.clear();
}

/** Execute the leap process */
int LeapInterface::execute_leap(FileName const& input) const {
# ifdef _MSC_VER
  return 1;
# else /* _MSC_VER */
  // Add source commands and quit to temp leap input
  CpptrajFile leapin;

  if (leapin.OpenWrite(input)) {
    mprinterr("Error: Could not open temporary leap input file '%s'\n", input.full());
    return 1;
  }

  for (Sarray::const_iterator it = input_files_.begin(); it != input_files_.end(); ++it)
    leapin.Printf("source %s\n", it->c_str());
  for (Sarray::const_iterator it = commands_.begin(); it != commands_.end(); ++it)
    leapin.Printf("%s\n", it->c_str());
  leapin.Printf("quit\n");

  leapin.CloseFile();

  std::string cmd("tleap -f " + input.Full());
  if (debug_ > 0) mprintf("DEBUG: %s\n", cmd.c_str());

//  int err = system( cmd.c_str() );
  Sarray errorMessages;

  CpptrajFile leapout;
  if (debug_ > 0 || !leapOutName_.empty()) {
    if (leapout.OpenWrite(leapOutName_)) {
      mprinterr("Error: Could not open leap output file '%s'\n", leapOutName_.c_str());
      return 1;
    }
  }

  FILE* file = popen(cmd.c_str(), "r");
  if (file == 0) {
    mprinterr("Error: Could not execute '%s'\n", cmd.c_str());
    return 1;
  }

  static const unsigned int BUFSIZE = 1023;
  char buffer[BUFSIZE+1];
  char* ptr = fgets(buffer, BUFSIZE, file);
  bool cleanExit = false;
  std::string exitLine;
  while (ptr != 0) {
    std::string line(ptr);
    if (leapout.IsOpen())
      leapout.Write(line.c_str(), line.size());
    std::size_t found = line.find("FATAL:");
    if (found != std::string::npos)
      errorMessages.push_back( NoTrailingWhitespace(line) );
    found = line.find("Exiting LEaP");
    if (found != std::string::npos) {
      cleanExit = true;
      exitLine = line;
    }
    ptr = fgets(buffer, BUFSIZE, file);
  }

  leapout.CloseFile();
  pclose(file);

  if (!cleanExit) return 1;
  if (debug_ > 0) mprintf("DEBUG: Leap Exit line '%s'\n", exitLine.c_str());
  ArgList exitArgs(exitLine, " :=;.");
  int nerr = exitArgs.getKeyInt("Errors", -1);
  int nwarn = exitArgs.getKeyInt("Warnings", -1);
  int nnotes = exitArgs.getKeyInt("Notes", -1);
  mprintf("\tLeap Errors= %i  Warnings= %i  Notes= %i\n", nerr, nwarn, nnotes);
  if (nerr == -1 || nerr > 0) {
    for (Sarray::const_iterator it = errorMessages.begin(); it != errorMessages.end(); ++it)
      mprintf("%s\n", it->c_str());
    return 1;
  }
  return 0;
# endif /* _MSC_VER */
}

/** Run leap. */
int LeapInterface::RunLeap() const {
# ifdef _MSC_VER
  // popen/pclose does not really work on windows.
  mprinterr("Error: LEaP interface cannot be used on windows.\n");
  return 1;
# else
  // Create temporary leap input file
  FileName tmp_leap_input = File::GenTempName();

  int err = execute_leap( tmp_leap_input );

  // Free up temp name
  File::FreeTempName( tmp_leap_input );
# endif
  return err;
}
