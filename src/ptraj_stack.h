#ifndef INC_PTRAJ_STACK_H
#define INC_PTRAJ_STACK_H
// ptraj - stackType
// originally from utility.h

typedef struct _stackType {
  void *entry;
  struct _stackType *next;
} stackType;


//typedef void (*fxnPrintStackEntry)(void *);

void pushBottomStack( stackType **, void * );
void pushStack( stackType **, void * );
void *popStack( stackType ** );
void clearStack( stackType ** );
//void printStack( stackType **, fxnPrintStackEntry, char *);
#endif
