// ptraj - stackType
// originally from utility.h

typedef struct _stackType {
  void *entry;
  struct _stackType *next;
} stackType;

void pushBottomStack( stackType **, void * );
void pushStack( stackType **, void * );
