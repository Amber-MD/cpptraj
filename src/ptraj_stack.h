#ifndef INC_PTRAJ_STACK_H
#define INC_PTRAJ_STACK_H
// ptraj - stackType
// originally from utility.h

typedef struct _stackType {
  void *entry;
  struct _stackType *next;
} stackType;


// GLOBAL STACK VARIABLES
#ifdef PTRAJ_STACK_MODULE
stackType *vectorStack = NULL;
stackType *matrixStack = NULL;
stackType *modesStack = NULL;
#else
extern stackType *vectorStack;
extern stackType *matrixStack;
extern stackType *modesStack;
#endif
//typedef void (*fxnPrintStackEntry)(void *);

void pushBottomStack( stackType **, void * );
void pushStack( stackType **, void * );
void *popStack( stackType ** );
void clearStack( stackType ** );
//void printStack( stackType **, fxnPrintStackEntry, char *);
#endif
