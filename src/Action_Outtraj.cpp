// Outtraj 
#include "Action_Outtraj.h"
#include "CpptrajStdio.h"
#include "PtrajMpi.h" //worldsize

// CONSTRUCTOR
Outtraj::Outtraj() {
  //fprintf(stderr,"Outtraj Con\n");
} 

// DESTRUCTOR
Outtraj::~Outtraj() { }

/*
 * Outtraj::init()
 * Action wrapper for trajout
 */
int Outtraj::init() {

  mprintf("    OUTTRAJ: [%s]\n",A->ArgLine());
  return ( outtraj.Add(A,PFL,worldsize) );
} 

/*
 * Outtraj::setup()
 * Output trajectory setup done right before first write
 */
//int Outtraj::setup() {
//  return 0;  
//}

/*
 * Outtraj::action()
 */
int Outtraj::action() {

  return ( outtraj.Write(currentFrame, F, P) );
} 


