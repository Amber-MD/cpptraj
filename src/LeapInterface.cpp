#include "LeapInterface.h"
#include "CpptrajFile.h"
#include "CpptrajStdio.h"
#include "File_TempName.h"
#include <cstdio>

using namespace Cpptraj;
/** CONSTRUCTOR */
LeapInterface::LeapInterface() {}

/** Set up leap to run. */
int LeapInterface::AddInputFile(std::string const& fname) {
  if (fname.empty()) return 1;
  input_files_.push_back( fname );
  return 0;
}

/** Clear all input files. */
void LeapInterface::ClearInputFiles() {
  input_files_.clear();
}

/** Execute the leap process */
int LeapInterface::execute_leap(FileName const& input) const {
  // Add source commands and quit to temp leap input
  CpptrajFile leapin;

  if (leapin.OpenWrite(input)) {
    mprinterr("Error: Could not open temporary leap input file '%s'\n", input.full());
    return 1;
  }

  for (Sarray::const_iterator it = input_files_.begin(); it != input_files_.end(); ++it)
    leapin.Printf("source %s\n", it->c_str());
  leapin.Printf("quit\n");

  leapin.CloseFile();

  CpptrajFile leapout;
  if (leapout.OpenWrite(leapOutName_)) {
    mprinterr("Error: Could not open leap output file '%s'\n", leapOutName_.c_str());
    return 1;
  }

  std::string cmd("tleap -f " + input.Full());
  FILE* file = popen(cmd.c_str(), "r");
  if (file == 0) {
    mprinterr("Error: Could not execute '%s'\n", cmd.c_str());
    return 1;
  }

  static const unsigned int BUFSIZE = 1023;
  char buffer[BUFSIZE+1];
  char* ptr = fgets(buffer, BUFSIZE, file);
  while (ptr != 0) {
    std::string line(ptr);
    leapout.Write(line.c_str(), line.size());
    ptr = fgets(buffer, BUFSIZE, file);
  }

  leapout.CloseFile();
  pclose(file);

  return 0;
}

/** Run leap. */
int LeapInterface::RunLeap() const {
  // Create temporary leap input file
  FileName tmp_leap_input = File::GenTempName();

  int err = execute_leap( tmp_leap_input );

  // Free up temp name
  File::FreeTempName( tmp_leap_input );
  return err;
}
