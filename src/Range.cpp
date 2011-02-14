// Range
#include <cstdlib>
#include <cstring>
#include "ArgList.h"
#include "CpptrajStdio.h"
#include "Range.h"

// CONSTRUCTOR
Range::Range() {
  rangeArg=NULL;
}

// DESTRUCTOR
Range::~Range() {
  if (rangeArg!=NULL) free(rangeArg);
}

/*
 * Range::SetRange()
 * Given an argument containing numbers separated by , and -, construct
 * a list of numbers corresponding to the argument.
 */
int Range::SetRange(char *ArgIn) {
  char *temp, *arg;
  int R[2], upper, err;
  ArgList *CommaList, *DashList;
  std::list<int>::iterator it;

  //mprintf("DEBUG: SetRange(%s)\n",ArgIn);

  if (ArgIn==NULL) return 1;
  this->clear();

  // Set rangeArg - realloc to allow re-use
  rangeArg=(char*) realloc( rangeArg, (strlen(ArgIn)+1) * sizeof(char));
  strcpy(rangeArg,ArgIn);
  // Copy ArgIn to temp location to avoid modifying it with strtok
  temp=(char*) malloc( (strlen(ArgIn)+1) * sizeof(char));
  strcpy(temp,ArgIn);
  // Split range by comma
  CommaList = new ArgList(temp, ",");
  CommaList->ResetAll();
  err=0;
  while ( (arg = CommaList->getNextString())!=NULL ) {
    // Then split by dash
    DashList = new ArgList(arg, "-");
    DashList->ResetAll();
    R[0] = DashList->getNextInteger(-1);
    R[1] = DashList->getNextInteger(-1);
    delete DashList;
    if (R[0]==-1) {
      mprintf("Error: Range::SetRange(%s): Range is -1 for %s\n",ArgIn,DashList->ArgLine());
      err=1;
      break;
    }
    upper = R[1];
    if (upper==-1) upper=R[0];
    upper++; // Want up to and including the upper argument
    if ( this->SetRange(R[0], upper) )
      mprintf("Warning: Converting %s to range [%i-%i] is not valid.\n",ArgIn,R[0],R[1]);
  }

  free(temp);
  delete CommaList;
  // Dont return an empty list
  if ( err>0 || this->empty() ) 
    return 1;
  
  // Sort frames using default comparison
  this->sort();
  //for (it=this->begin(); it!=this->end(); it++)
  //  fprintf(stdout,"RangeList= %i\n",*it); 
  // Remove duplicates
  err=-1;
  //fprintf(stdout,"Size of RangeList is %lu\n",this->size());
  it=this->begin();
  while (!this->empty()) {
    //fprintf(stdout,"     List= %i  Last= %i",*it,err);
    upper=*it;
    if (*it == err) {
      //fprintf(stdout,"REMOVING."); 
      // Erasing effectively increments the iterator
      it=this->erase(it);
    } else {
      it++;
    }
    err=upper;
    //fprintf(stdout,"\n");
    // If we are past the last element exit now
    if (it == this->end()) break;
  }

  return 0;
}

/*
 * Range::SetRange()
 * Given a start and end number, set up a range from start to (not 
 * including) end.
 * Check that end is greater than start so that the range list is
 * not empty.
 */
int Range::SetRange(int start, int end) {
  int range;

  if (end <= start) {
    mprintf("Error: Range::SetRange: end (%i) <= start (%i)\n",end,start);
    return 1;
  }

  for (range=start; range < end; range++)
    this->push_back(range);

  return 0;
}

/*
 * Range::RangeArg()
 * Return the range argument
 */
char *Range::RangeArg() {
  return rangeArg;
}

/*
 * Range::PrintRange()
 * Print all numbers in the range to a line. Increment by offset.
 */
void Range::PrintRange(const char* header, int offset) {
  if (header!=NULL)
    mprintf("%s",header);
  for (std::list<int>::iterator it=this->begin(); it!=this->end(); it++)
    mprintf(" %i",(*it)+offset);
  mprintf("\n");
}
 
