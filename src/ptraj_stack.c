// ptraj_stack
#include "ptraj_stack.h"
#include <stdlib.h>

// pushStack
void pushStack( stackType **stackp, void *entry ) {
  stackType *sp;

  sp = (stackType*) malloc( sizeof(stackType) );
  sp->entry = entry;
  sp->next = *stackp;
  *stackp = sp;
}   
  
// pushBottomStack  
void pushBottomStack( stackType **stackp, void *entry ) {   
  stackType *sp;
  
  if ( *stackp == NULL )
    pushStack( stackp, entry );
  else {
    sp = *stackp;
    while ( sp->next != NULL ) sp = sp->next;
    sp->next = (stackType*) malloc( sizeof(stackType) );
    sp->next->entry = entry;
    sp->next->next = NULL;
  }
}


