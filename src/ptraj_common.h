#ifndef INC_PTRAJ_COMMON_H
#define INC_PTRAJ_COMMON_H
// ptraj_common.h
#ifndef BUFFER_SIZE
#define BUFFER_SIZE 1024
#endif

// GLOBAL DEBUG LEVEL
#ifdef PTRAJ_COMMON_MODULE
int prnlev = 0;
#else
extern int prnlev;
#endif

void error(char *, char *, ...);
void warning(char *, char *, ...);
void *safe_malloc(size_t);
void *safe_realloc(void *, size_t , size_t);
void safe_free(void *);
FILE *safe_fopen(char *, char *);
void safe_fclose(FILE *);

#endif
