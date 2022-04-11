#include <cctype> // isspace
#include "CmdInput.h"
#include "StringRoutines.h" // RemoveTrailingWhitespace()

/** '#' indicates the beginning of a comment, backslash at the end of a line
  * indicates continuation.
  * \return 1 if more lines need to be read for this input.
  * \return 0 if no more input needed.
  * \return -1 if an error occurred.
  */
int CmdInput::AddInput(const char* lineIn) {
  if (lineIn == 0 || *lineIn == '\0') return 0;
  // Skip leading whitespace
  const char* ptr = lineIn;
  while ( isspace(*ptr) && *ptr != '\0' )
    ++ptr;
  if (*ptr == '\0') return 0;
  std::string line( ptr );
  // Remove trailing whitespace
  RemoveTrailingWhitespace( line );
  if (line.empty()) return 0;
  // Terminal backslash requests a continuation of input.
  size_t end = line.size();
  bool needMoreInput = (line[end-1] == '\\');
  if (needMoreInput) {
    size_t bs_pos = line.find_last_of('\\');
    line.erase( bs_pos, 1 );
  }
  // Find '#' not preceded by blackslash; indicates comment.
  // Remove it and all after.
  end = line.find_first_of('#');
  if (end != std::string::npos) {
    if (end == 0 || line[end-1] != '\\')
      line.erase( line.begin() + end, line.end() );
  }
  if (line.empty()) {
    if (needMoreInput) return 1;
    return 0;
  }
  // Add to current input, removing any consecutive whitespace.
  // TODO Convert all tabs to space?
  input_ += line[0];
  unsigned int idx1 = 1;
  unsigned int idx0 = input_.size() - 1;
  for (; idx1 != line.size(); ++idx1)
    if (!isspace( line[idx1] )) {
      input_ += line[idx1];
      idx0++;
    } else if (!isspace( input_[idx0] )) {
      input_ += line[idx1];
      idx0++;
    }
  if (needMoreInput) return 1;
  return 0;
}
