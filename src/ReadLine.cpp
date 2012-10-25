#include <cstdio> // required for readline
#include <cstdlib> // free
#define READLINE_LIBRARY
#include <readline.h>
#include <history.h>
#include "ReadLine.h"

void ReadLine::GetInput() {
  char* line = readline("> ");
  if (line && *line) add_history(line);
  input_.assign( line );
  free(line);
}
  
