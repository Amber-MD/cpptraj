// ptraj_stack
#include "ptraj_stack.h"
#include <stdlib.h>

// pushStack()
void pushStack( stackType **stackp, void *entry ) {
  stackType *sp;

  sp = (stackType*) malloc( sizeof(stackType) );
  sp->entry = entry;
  sp->next = *stackp;
  *stackp = sp;
}   
  
// pushBottomStack()  
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

// popStack()
void *popStack( stackType **stackp ) {
  void *entry;
  stackType *sp;

  if ( *stackp == NULL ) {
    return( (char *) NULL );
  }

  sp = *stackp;
  entry = sp->entry;

  *stackp = sp->next;
  sp->next=NULL;
  free( sp );
  return( entry );
}

// clearStack()
void clearStack( stackType **stackp ) {
  stackType *sp, *tmp;

  if (stackp == NULL) {
    return;
  }

  tmp = NULL;
  for (sp = *stackp; sp != NULL; ) {
    if ( tmp != NULL ) free( (void *) tmp);
    free( (void *) sp->entry);
    sp->entry = NULL;
    tmp = sp;
    sp = sp->next;
  }
  *stackp = NULL;
}
