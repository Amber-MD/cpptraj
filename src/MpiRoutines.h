#ifndef INC_MPIROUTINES_H
#define INC_MPIROUTINES_H
#ifdef __cplusplus
extern "C" {
#endif
/*! \file MpiRoutines.h
    \brief Cpptraj interface to C MPI routines.

    MPIROUTINES_MODULE is defined in MpiRoutines.c. For all other files
    including MpiRoutines.h worldrank and worldsize should be extern.
 */
// If debug is not defined stdio.h will be included in MpiRoutines.c
#ifdef DEBUG
#  include <stdio.h>
#endif
#include <sys/types.h> // off_t

#ifdef MPIROUTINES_MODULE
int worldrank;
int worldsize;
#  ifdef DEBUG
FILE *mpidebugfile;
#  endif
#else
extern int worldrank;
extern int worldsize;
#  ifdef DEBUG
extern FILE *mpidebugfile;
#  endif
#endif

/// This allows abstraction of the MPI_File type so no other files need mpi.h
typedef struct parallelStructType *parallelType;

#ifdef MPI
/* ========== Routines that require MPI ========== */
void printMPIerr(int, const char *);
int checkMPIerr(int, const char *);
int parallel_check_error(int );
#endif

/* ========== Routines that do not require MPI ========== */
void dbgprintf(const char *, ...);
int parallel_init(int, char **);
int parallel_end();
void parallel_barrier();
int parallel_sum(double *, double *, int );
int parallel_openFile_read(parallelType, const char *);
int parallel_open_file_write(parallelType, const char *);
int parallel_closeFile(parallelType);
int parallel_fread(parallelType, void*, int );
int parallel_fwrite(parallelType, void *, int);
int parallel_fseek(parallelType, off_t, int);
char *parallel_fgets(parallelType, char *, int);
int parallel_setSize(parallelType, long int);
int parallel_sendMaster(void *, int, int, int);
int parallel_allreduce(void *, void *, int, int, int);

#ifdef __cplusplus
}
#endif
#endif
