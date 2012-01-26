// Range
#include <cstring>
#include "Range.h"
#include "ArgList.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Range::Range() { }

// DESTRUCTOR
Range::~Range() { }

// COPY CONSTRUCTOR
Range::Range(const Range &rhs) {
  rangeArg = rhs.rangeArg;
  rangeList = rhs.rangeList;
  // Dont copy over the iterator, set to beginning
  rangeNum = rangeList.begin();
}

// ASSIGNMENT OPERATOR
Range &Range::operator=(const Range &rhs) {
  // Check for self assignment
  if ( this == &rhs ) return *this;
  rangeArg = rhs.rangeArg;
  rangeList = rhs.rangeList;
  // Dont copy over the iterator, set to beginning
  rangeNum = rangeList.begin();
  // Return *this
  return *this;
}

// Range::SetRange()
/** Given an argument containing numbers separated by "," (concatentation), and 
  * "-" (number range), construct an ordered list of numbers corresponding to 
  * the argument. Remove any duplicate numbers.
  * \return 0 on success, 1 on error.
  */
int Range::SetRange(char *ArgIn) {
  char *temp, *arg;
  int R[2], upper, err;
  ArgList CommaList, DashList;
  std::list<int>::iterator it;

  //mprintf("DEBUG: SetRange(%s)\n",ArgIn);

  if (ArgIn==NULL) return 1;
  rangeList.clear();

  // Set rangeArg
  rangeArg.assign(ArgIn);
  // Check if ArgIn is a mask expression
  if ( strchr( ArgIn, ':')!=NULL ||
       strchr( ArgIn, '@')!=NULL ||
       strchr( ArgIn, '*')!=NULL   )
  {
    mprinterr("Error: Using a mask expression for range (%s)\n",ArgIn);
    mprinterr("       Ranges only contain digits, dashes, and commas (e.g. 3-5,8-10)\n");
    return 1;
  } 
  // Copy ArgIn to temp location to avoid modifying it with strtok
  temp = new char[ strlen(ArgIn)+1 ];
  strcpy(temp,ArgIn);
  // Split range by comma
  CommaList.SetList(temp, ",");
  err=0;
  while ( (arg = CommaList.getNextString())!=NULL ) {
    // Then split by dash
    DashList.SetList(arg, "-");
    R[0] = DashList.getNextInteger(-1);
    R[1] = DashList.getNextInteger(-1);
    if (R[0]==-1) {
      mprintf("Error: Range::SetRange(%s): Range is -1 for %s\n",ArgIn,DashList.ArgLine());
      err=1;
      break;
    }
    upper = R[1];
    if (upper==-1) upper=R[0];
    upper++; // Want up to and including the upper argument
    if ( this->SetRange(R[0], upper) )
      mprintf("Warning: Converting %s to range [%i-%i] is not valid.\n",ArgIn,R[0],R[1]);
  }

  delete[] temp;
  // Dont return an empty list
  if ( err>0 || rangeList.empty() ) 
    return 1;
  
  // Sort frames using default comparison
  rangeList.sort();
  //for (it=rangeList.begin(); it!=rangeList.end(); it++)
  //  fprintf(stdout,"RangeList= %i\n",*it); 
  // Remove duplicates
  err=-1;
  //fprintf(stdout,"Size of RangeList is %lu\n",rangeList.size());
  it=rangeList.begin();
  while (!rangeList.empty()) {
    //fprintf(stdout,"     List= %i  Last= %i",*it,err);
    upper=*it;
    if (*it == err) {
      //fprintf(stdout,"REMOVING."); 
      // Erasing effectively increments the iterator
      it=rangeList.erase(it);
    } else {
      it++;
    }
    err=upper;
    //fprintf(stdout,"\n");
    // If we are past the last element exit now
    if (it == rangeList.end()) break;
  }

  // Set iterator to beginning
  this->Begin();

  return 0;
}

// Range::SetRange()
/** Given a start and end number, set up a range from start to (not 
  * including) end.  
  */
int Range::SetRange(int start, int end) {
  int range;
  
  //Check that end is greater than start so that the range list is
  if (end <= start) {
    mprintf("Error: Range::SetRange: end (%i) <= start (%i)\n",end,start);
    return 1;
  }

  for (range=start; range < end; range++)
    rangeList.push_back(range);

  // Set iterator to beginning
  this->Begin();

  return 0;
}

// Range::SetRange()
/** Assign this range from another range. */
// NOTE: Obsolete?
void Range::SetRange(Range *rhs) {
  int num;

  rhs->Begin();
  while (rhs->NextInRange(&num))
    this->AddToRange(num);

  // Set iterator to beginning
  this->Begin();
}

// Range::ShiftBy()
/** Shift all numbers in range by specified value. */
void Range::ShiftBy(int val) {
  for (rangeNum = rangeList.begin(); rangeNum != rangeList.end(); rangeNum++)
    *rangeNum += val;
  // Set iterator to beginning
  this->Begin();  
}

// Range::AddToRange()
/** Add a number to the range. Range is NOT explicitly sorted in this case.
  */
void Range::AddToRange(int num) {
  rangeList.push_back(num);
}

// Range::Begin()
/** Set iterator to the beginning of the list.
  */
void Range::Begin() {
  rangeNum = rangeList.begin();
}

// Range::NextInRange()
/** Set num to the current number in the range and increment the iterator. 
  * Return true if number was set, false if no more numbers in range.
  */
bool Range::NextInRange(int *num) {
  // NOTE: No check for NULL here.
  if (rangeNum == rangeList.end()) return false;
  *num = *rangeNum;
  rangeNum++;
  return true;
}

// Range::End()
/** \return true if at the end of the range.
  */
bool Range::End() {
  if (rangeNum == rangeList.end()) return true;
  return false;
}

// Range::Next()
/** Increment the iterator. */
void Range::Next() {
  rangeNum++;
}

// Range::Current()
/** Return the current number in the range. */
int Range::Current() {
  return *rangeNum;
}

// Range::RemoveFromRange()
/** Remove all instances of num from the range. */
void Range::RemoveFromRange(int num) {
  std::list<int>::iterator it=rangeList.begin();
  while (it!=rangeList.end()) {
    if (*it == num) 
      it = rangeList.erase(it);
    else
      it++;
  }
}

// Range::RangeArg()
/** Return the range argument */
char *Range::RangeArg() {
  return (char*)rangeArg.c_str();
}

// Range::PrintRange()
/** Print all numbers in the range to a line. Increment by offset. */
void Range::PrintRange(const char* header, int offset) {
  if (header!=NULL)
    mprintf("%s",header);
  for (std::list<int>::iterator it=rangeList.begin(); it!=rangeList.end(); it++)
    mprintf(" %i",(*it)+offset);
  //mprintf("\n");
}
 
