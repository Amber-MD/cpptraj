// Outtraj 
#include "Action_Outtraj.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Outtraj::Outtraj() {
  //fprintf(stderr,"Outtraj Con\n");
  min=0.0;
  max=0.0;
  Dset=NULL;
} 

// DESTRUCTOR
Outtraj::~Outtraj() { }

/*
 * Outtraj::init()
 * Action wrapper for trajout.
 * Expected call: outtraj <filename> [ trajout args ] 
 *                        [maxmin <dataset> min <min> max <max>
 */
int Outtraj::init() {
  char *datasetName;

#ifdef MPI
  mprintf("ERROR: OUTTRAJ currently not functional with MPI.\n");
  return 1;
#endif

  mprintf("    OUTTRAJ: Will write to [%s]\n",A->Arg(1));
  outtraj.SetDebug(debug);
  // If maxmin, get the name of the dataset as well as the max and min values.
  datasetName = A->getKeyString("maxmin",NULL);
  if (datasetName!=NULL) {
    Dset = DSL->Get(datasetName);
    if (Dset==NULL) {
      mprintf("Error: Outtraj maxmin: Could not get dataset %s\n",datasetName);
      return 1;
    } else {
      // Currently dont allow for string datasets
      if (Dset->Type()==STRING) {
        mprintf("Error: Outtraj maxmin: String dataset (%s) not supported.\n",datasetName);
        return 1;
      }
      max = A->getKeyDouble("max",0.0);
      min = A->getKeyDouble("min",0.0);
      mprintf("             maxmin: Printing trajectory frames based on %lf <= %s <= %lf\n",
              min, datasetName, max);
    }
  }

  return ( outtraj.Add(A,PFL) );
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
  double dVal;
  int iVal;

  // If dataset defined, check if frame is within max/min
  if (Dset!=NULL) {
    if (Dset->Type() == DOUBLE) {
      if (Dset->Get(&dVal, currentFrame)) return 1;
    } else if (Dset->Type() == INT) {
      if (Dset->Get(&iVal, currentFrame)) return 1;
      dVal = (double) iVal;
    } else
      return 1;
    //mprintf("DBG: maxmin: dVal = %lf\n",dVal);
    // If value from dataset not within min/max, exit now.
    if (dVal < min || dVal > max) return 0;
  }
  if ( outtraj.Write(currentFrame, F, P) != 0 ) return 1;
  return 0;
} 

/*
 * Outtraj::print()
 * Close trajectory.
 */
void Outtraj::print() {
  mprintf("  OUTTRAJ: [%s] Wrote %i frames.\n",A->Arg(1),outtraj.front()->CurrentFrame());
  outtraj.Close();
}

