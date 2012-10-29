#include <cstdio> // required for readline
#include <cstdlib> // free
#include <cstring>
#define READLINE_LIBRARY
#include <readline.h>
#include <history.h>
#include "ReadLine.h"

/** Get next input line with readline. Lines terminated with a backslash
  * will be concatenated. Comments will be ignored.
  */
int ReadLine::GetInput() {
  input_.clear();
  char* line = readline("> ");
  if (line == 0) return 1; // EOF
  input_ += line;
  // Terminal backslash requests a continuation of the line
  size_t end = strlen( line );
  while (end > 1 && line[end - 1] == '\\') {
    // Remove that backlash
    size_t bs_pos = input_.find_last_of('\\');
    input_.erase( bs_pos, 1 );
    free( line );
    line = readline("");
    if (line == 0) break;
    input_ += line;
    end = strlen( line );
  }
  // Remove leading whitespace.
  std::string::iterator beg = input_.begin();
  while ( beg != input_.end() && isspace(*beg) )
    beg = input_.erase(beg);
  // Find '#' not preceded by blackslash; indicates comment.
  // Remove it and all after.
  end = input_.find_first_of('#');
  if (end != std::string::npos) {
    if (end == 0 || (end > 0 && input_[end-1] != '\\'))
      input_.erase( input_.begin() + end, input_.end() );
  }
  // Add line to history
  if (!input_.empty()) add_history(input_.c_str());
  if (line != 0) free( line );
  return 0;
}
