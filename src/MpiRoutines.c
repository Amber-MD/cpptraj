/*! \file MpiRoutines.c
    \brief Cpptraj interface to MPI C routines.
 */
#define MPIROUTINES_MODULE
/* MPIROUTINES_MODULE is defined here so that worldrank and worldsize are defined
 * as int in the MpiRoutines header file.
 */
#include "MpiRoutines.h"
// If DEBUG is defined stdio.h will be included in MpiRoutines.h
#ifndef DEBUG
#  include <stdio.h>
#endif
#include <stdlib.h>
#ifdef MPI
#  include "mpi.h"
#endif

/// This allows abstraction of the MPI_File type so no other files need mpi.h
struct parallelStructType {
#ifdef MPI
  MPI_File* mfp;
#else
  void* mfp;
#endif
};

#ifdef MPI
/* ========== Routines that require MPI ========== */
// printMPIerr()
/** Wrapper for MPI_Error string.  */
void printMPIerr(int err, const char *routineName) {
  int len,eclass,i;
  char buffer[1024];

  MPI_Error_string(err,buffer,&len);
  MPI_Error_class(err,&eclass);
  // Remove newlines from MPI error string
  for (i=0; i<len; i++)
    if (buffer[i]=='\n') buffer[i]=':';
  fprintf(stdout,"[%i] MPI ERROR %d: %s: [%s]\n",worldrank,eclass,routineName,buffer);

  return;
}

// checkMPIerr()
/** \return 1 and print error message if MPI action failed. */
int checkMPIerr(int err, const char *routineName) {
  if (err!=MPI_SUCCESS) {
    printMPIerr(err, routineName);
    return 1;
  } 
  return 0;
}

// parallel_check_error()
/** All ranks pass in error value. Compute the sum. return non-zero if error on
  * any nodes.
  */
int parallel_check_error(int err) {
  int errtotal;

  errtotal=0;
  MPI_Allreduce(&err,&errtotal,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  return errtotal;
}
#endif

/* ========== Routines that do not require MPI ========== */
#ifdef DEBUG
// dbgprintf()
/** Print to mpidebugfile */
void dbgprintf(const char *format, ...) {
  va_list args;

  va_start(args,format);
  //fprintf(mpidebugfile,"[%i] ",worldrank);
  vfprintf(mpidebugfile,format,args);
  fflush(mpidebugfile);
  va_end(args);
  return;
}

// parallel_debug_init()
/** Open a file named Thread.worldrank for this thread */
int parallel_debug_init() {
  char outfilename[32];

  // DEBUG
  sprintf(outfilename,"Thread.%03i",worldrank);
  mpidebugfile=fopen(outfilename,"w");
  if (mpidebugfile==NULL) {
    fprintf(stderr,"[%i]\tCould not open debug file:\n",worldrank);
    perror("");
    return 1;
  } /*else { 
    dbgprintf("MPI DEBUG: %s %p\n",outfilename,mpidebugfile);
    fprintf(stdout,"MPI DEBUG: %s %p\n",outfilename,mpidebugfile);
  }*/
  return 0;
}

// parallel_debug_end()
/** Close Thread.worldrank file.  */
int parallel_debug_end() {
  if (mpidebugfile!=NULL)
    fclose(mpidebugfile);
  return 0;
}

#endif

// parallel_init()
/** Initialize MPI */
int parallel_init(int argc, char **argv) {
#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &worldsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldrank);
#else
  worldrank=0;
  worldsize=1;
#endif
#ifdef DEBUG
  parallel_debug_init();
#endif
  return 0;
}

// parallel_end()
/** Finalize MPI */
int parallel_end() {
#ifdef DEBUG
  parallel_debug_end();
#endif
#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif
  return 0;
}

// parallel_barrier()
/** Use MPI barrier.  */
void parallel_barrier() {
#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

// parallel_sum()
/** Use MPI_REDUCE to sum up the values in sendbuffer and place them in
  * recvbuffer on master.
  */
int parallel_sum(double *recvBuffer, double *sendBuffer, int N) {
#ifdef MPI
  int err;

  err=MPI_Reduce(sendBuffer, recvBuffer, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if ( parallel_check_error(err) ) {
    checkMPIerr(err,"parallel_sum");
    return 1;
  }
#endif
  return 0;
}

// parallel_openFile_read()
/** Use MPI to open a file for reading.  */
int parallel_openFile_read(parallelType pfile, const char *filename) {
#ifdef MPI
  int err;
  MPI_File *mfp;

  pfile->mfp=NULL;
  //if (prnlev>0) 
    fprintf(stdout,"[%i]\tparallel_openFile_read: Opening input file %s\n",worldrank,filename);
  mfp = (MPI_File*) malloc(sizeof(MPI_File));
  err=MPI_File_open(MPI_COMM_WORLD, (char*) filename, MPI_MODE_RDONLY, MPI_INFO_NULL, mfp);
  if (err!=MPI_SUCCESS) printMPIerr(err,"parallel_openFile_read()");
  // Check that everyone opened the file. 
  if (parallel_check_error(err)!=0) {
    free(mfp);
    return 1;
  } else {
    pfile->mfp=mfp;
    return 0;
  }
#endif
  return 1;
}

// parallel_open_file_write()
/** Use MPI to open a file for writing, delete if present. */
int parallel_open_file_write(parallelType pfile, const char *filename) {
#ifdef MPI
  int err;
  MPI_File *mfp;
  
  pfile->mfp=NULL;
  // Remove file if present
  MPI_File_delete((char*)filename,MPI_INFO_NULL);

#  ifdef DEBUG 
    dbgprintf("\tparallel_open_file_write: Opening output file %s\n",filename);
#  endif
  mfp = (MPI_File*) malloc(sizeof(MPI_File));
  err=MPI_File_open(MPI_COMM_WORLD, (char*) filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, 
                    MPI_INFO_NULL, mfp);
  if (err!=MPI_SUCCESS) printMPIerr(err,"parallel_open_file_write()");
  // Check that everyone opened the file. 
  if (parallel_check_error(err)!=0) {
    free(mfp);
    return 1;
  } else {
    pfile->mfp=mfp;
    return 0;
  }
#endif
  return 1;
}

// parallel_closeFile()
/** Close MPI file.  */
int parallel_closeFile(parallelType pfile) {
#ifdef MPI
  int err;
  if (pfile->mfp==NULL) return 0;
  //MPI_Barrier(MPI_COMM_WORLD);
#  ifdef DEBUG
    dbgprintf("\tparallel_closeFile: Closing file.\n");
#  endif
  err=MPI_File_close(pfile->mfp);
  if (err!=MPI_SUCCESS) printMPIerr(err,"parallel_closeFile");
  // Free the memory used by the file pointer
  free(pfile->mfp);
  pfile->mfp=NULL;
#endif 
 return 0;
}

// parallel_fread()    
/** fread using MPI routines.  */
int parallel_fread(parallelType pfile, void *buffer, int count) {
#ifdef MPI
  int err, actualCount;
  MPI_Status status;

  err=MPI_File_read( *(pfile->mfp),buffer,count,MPI_CHAR,&status);
  if (err!=MPI_SUCCESS) {
    printMPIerr(err,"parallel_fread");
    return -1;
  }
  err=MPI_Get_count(&status, MPI_CHAR, &actualCount);
  if (err!=MPI_SUCCESS) {
    printMPIerr(err,"parallel_fread");
    return -1;
  }
  if (actualCount==MPI_UNDEFINED) {
    fprintf(stdout,"[%i] Error in parallel_fread: Number of bytes read was undefined.\n",worldrank);
    return -1; 
  }

  return actualCount;
#endif
  return -1;
}

// parallel_fwrite()
/** fwrite using MPI routines.  */
int parallel_fwrite(parallelType pfile, void *buffer, int count) {
#ifdef MPI
  int err;
  MPI_Status status;

#  ifdef DEBUG
  char *temp;
  temp=(char*) buffer;
  //dbgprintf("Calling MPI write(%i): [%s]\n",count,temp);
  dbgprintf("Calling MPI write(%i):\n",count);
#  endif

  err=MPI_File_write( *(pfile->mfp),buffer,count,MPI_CHAR,&status);
  if (err!=MPI_SUCCESS) {
    printMPIerr(err,"parallel_fwrite");
    return 1;
  }

  return 0;
#endif
  return 1;
}

// parallel_fseek()
int parallel_fseek(parallelType pfile, off_t offset, int origin) {
#ifdef MPI
  int err, org;
  MPI_Offset mpiOffset;
 
  org = MPI_SEEK_SET; 
  switch ( origin ) {
    case SEEK_SET : org=MPI_SEEK_SET; break;
    case SEEK_CUR : org=MPI_SEEK_CUR; break;
    case SEEK_END : org=MPI_SEEK_END; break;
    default: return 1;
  }

  mpiOffset = (MPI_Offset) offset;
 
  err=MPI_File_seek( *(pfile->mfp), mpiOffset, org);
  if (err!=MPI_SUCCESS) {
    printMPIerr(err,"parallel_fseek");
    return 1;
  } else 
    return 0;
#endif
  return 1;
}

// parallel_fgets()
/** Like fgets, use mpi file routines to get all chars up to and including 
  * null or newline. Returns buffer, or NULL on error. 
  */
char *parallel_fgets(parallelType pfile, char *buffer, int num) {
#ifdef MPI
  int i,err;
  
  for (i=0; i<num-1; i++) {
    err=MPI_File_read(*(pfile->mfp),buffer+i,1,MPI_CHAR,MPI_STATUS_IGNORE);
    if (err!=MPI_SUCCESS) {
      printMPIerr(err,"parallel_fgets");
      return NULL;
    }
    if (buffer[i]=='\n' || buffer[i]=='\0') {i++; break;} // Always have i be 1 more char ahead   
  } 
  // NOTE: Uncommenting the next 3 lines replaces newlines with NULL
  //if (i==num && buffer[i-1]=='\n')
  //  buffer[i-1]='\0';
  //else
    buffer[i]='\0';
    
  return buffer;
#endif
  return NULL;
} 

// parallel_setSize()
/** Set the size of mpi file to offset. */
int parallel_setSize(parallelType pfile, long int offset) {
#ifdef MPI
  int err;
  MPI_Offset mpiOffset;

  mpiOffset = (MPI_Offset) offset;

  err = MPI_File_set_size( *(pfile->mfp), mpiOffset);
  if (err!=MPI_SUCCESS) printMPIerr(err,"parallel_setSize");
  if (parallel_check_error(err)!=0) 
    return 1;

  return 0;
#endif
  return 1;
} 

// parallel_sendMaster()
/** If master    : receive specified values from rank.
  * If not master: send specified values from rank to master.
  * datatype: 0=MPI_INT
  *           1=MPI_DOUBLE
  */
int parallel_sendMaster(void *Buffer, int Count, int rank, int datatype) {
#ifdef MPI
  int err;
  MPI_Datatype currentType;

  if (worldsize==1) return 0;

  switch (datatype) {
    case 0 : currentType = MPI_INT; break;
    case 1 : currentType = MPI_DOUBLE; break;
    case 2 : currentType = MPI_CHAR; break;
    default: return 1;
  }

  if (worldrank>0) {
    // Non-master
    err = MPI_Send(Buffer, Count, currentType, 0, 1234, MPI_COMM_WORLD); 
    if (err!=MPI_SUCCESS) printMPIerr(err,"Sending data to master.");

  } else {
    // Master
    err = MPI_Recv(Buffer, Count, currentType, rank, 1234, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (err!=MPI_SUCCESS) printMPIerr(err,"Receiving data from non-master.");
  }

  if (parallel_check_error(err)!=0) return 1;

  return 0;
#endif
  return 1;
}

// parallel_allreduce()
/** Perform an mpi allreduce. */
int parallel_allreduce(void *Return, void *input, int count, int datatype, int op) {
#ifdef MPI
  int err;
  MPI_Datatype currentType;
  MPI_Op currentOp;

  if (worldsize==1) return 0;

  switch (datatype) {
    case 0 : currentType = MPI_INT; break;
    case 1 : currentType = MPI_DOUBLE; break;
    case 2 : currentType = MPI_CHAR; break;
    default: return 1;
  }  

  switch (op) {
    case 0 : currentOp = MPI_SUM; break;
    default: return 1;
  }

  err = MPI_Allreduce(input, Return, count, currentType, currentOp, MPI_COMM_WORLD);
  if (err!=MPI_SUCCESS) {
    printMPIerr(err, "Performing allreduce.\n");
    printf("[%i]\tError: allreduce failed for %i elements.\n",worldrank,count);
  }

  if (parallel_check_error(err)!=0) return 1;
  return 0;
#endif
  return 1;
}

