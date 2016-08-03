#include "CpptrajStdio.h"
#include <cstdio> // required for readline
#include <cstdlib> // free
#include <cstring>
#if defined (NO_READLINE) || defined (LIBCPPTRAJ)
#   include <stdexcept>
#   include <string>
#   include <iostream>
#else
#   define READLINE_LIBRARY
#   include <readline.h>
#   include <history.h>
#endif
#include "ReadLine.h"
#include "Command.h"

#if defined (NO_READLINE) || defined (LIBCPPTRAJ)
// CONSTRUCTOR
ReadLine::ReadLine() {}
#else

// duplicate_string()
static char* duplicate_string(const char* s) {
  char* r = (char*)malloc(strlen(s) + 1);
  strcpy(r, s);
  return r;
}

// command_generator()
/** Generator function for command completion.  STATE lets us know whether
  * to start from scratch; without any state (i.e. STATE == 0), then we
  * start at the top of the list.
  */
static char* command_generator(const char* text, int state) {
  static int list_index, len;
  const char *name;

  // If this is a new word to complete, initialize now. This includes
  // saving the length of TEXT for efficiency, and initializing the index
  // variable to 0.
  if (!state) {
    list_index = 0;
    len = strlen(text);
  }

  // Return the next name which partially matches from the command list.
  while ( (name = Command::CmdToken(list_index)) != 0 )
  {
    list_index++;
    if (strncmp(name, text, len) == 0)
      return (duplicate_string(name));
  }

  // If no names matched, then return NULL.
  return 0;
}

// cpptraj_completion()
static char** cpptraj_completion(const char* text, int start, int end) {
  char** matches = 0;
  // If this word is at the start of the line, assume it is a command.
  if (start == 0 || (strncmp(rl_line_buffer, "help ", 5)==0))
    matches = rl_completion_matches(text, command_generator);
  return matches;
}

// CONSTRUCTOR
ReadLine::ReadLine() {
  // Tell the completer that we want a crack first.
  rl_attempted_completion_function = cpptraj_completion;
}
#endif

// -----------------------------------------------------------------------------
/** Get next input line with readline. Lines terminated with a backslash
  * will be concatenated. Comments will be ignored.
  */
int ReadLine::GetInput() {
  input_.Clear();
# if defined (NO_READLINE) || defined (LIBCPPTRAJ)
  mprintf("> ");
  std::string inp;
  getline(std::cin, inp);
  if (inp.size() == 0 && std::cin.eof()) return 1; // EOF
  bool moreInput = input_.AddInput( inp.c_str() );
  while ( moreInput ) {
      getline(std::cin, inp);
      moreInput = input_.AddInput( inp.c_str() );
  }
# else
  char* line = readline("> ");
  if (line == 0) return 1; // EOF
  bool moreInput = input_.AddInput( line );
  while ( moreInput ) {
    free( line );
    line = readline("");
    moreInput = input_.AddInput( line );
  }
  // Add line to history
  if (!input_.Empty()) AddHistory(input_.str());
  if (line != 0) free( line );
# endif
  return 0;
}

void ReadLine::AddHistory(const char* line) {
# if defined (NO_READLINE) || defined (LIBCPPTRAJ)
  throw std::runtime_error("Internal Error: readline not enabled; history "
                           "storage not supported");
# else
  if (line != 0) add_history( line );
# endif
}

bool ReadLine::YesNoPrompt(const char* prompt) {
# if defined (NO_READLINE) || defined (LIBCPPTRAJ)
  mprintf(prompt);
  std::string line;
  getline(std::cin, line);
  if (line.size() == 0) return false;
# else
  char* line = readline(prompt);
  if (line == 0 || strlen( line ) < 1) return false;
# endif
  if (line[0] == 'y' || line[0] == 'Y') return true;
  return false;
}
