#include <cstdio> // required for readline
#include <cstdlib> // free
#define READLINE_LIBRARY
#include <readline.h>
#include <history.h>
#include "ReadLine.h"

int ReadLine::GetInput() {
  char* line = readline("> ");
  if (line) {
    if (*line) add_history(line);
  } else
    return 1;
  input_.assign( line );
  free(line);
  return 0;
}
  
