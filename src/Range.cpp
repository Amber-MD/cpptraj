// Range
#include <algorithm> // find
#include "Range.h"
#include "ArgList.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // ArrayToRangeExpression

/** CONSTRUCTOR */
Range::Range() { }

/** CONSTRUCTOR - Takes argument string as input. */
Range::Range( std::string const& argIn ) {
  if (!argIn.empty())
    SetRange( argIn );
}

/** CONSTRUCTOR - Single number. */
Range::Range(int start) { SetRange(start, start+1); }

/** CONSTRUCTOR - Range expression and offset. */
Range::Range( std::string const& argIn, int offsetIn) {
  if (!argIn.empty()) {
    SetRange( argIn );
    ShiftBy( offsetIn );
  }
}

/** COPY CONSTRUCTOR */
Range::Range(const Range &rhs) :
  rangeArg_(rhs.rangeArg_),
  rangeList_(rhs.rangeList_)
{}

/** ASSIGNMENT OPERATOR */
Range &Range::operator=(const Range &rhs) {
  // Check for self assignment
  if ( this == &rhs ) return *this;
  rangeArg_ = rhs.rangeArg_;
  rangeList_ = rhs.rangeList_;
  return *this;
}

// Range::SetRange()
/** Given an argument containing numbers separated by "," (concatentation), and 
  * "-" (number range), construct an ordered list of numbers corresponding to 
  * the argument. Remove any duplicate numbers.
  * \return 0 on success, 1 on error.
  */
int Range::SetRange(std::string const& ArgIn) {
  std::string arg;
  int R[2], upper;
  ArgList DashList;

  //mprintf("DEBUG: SetRange(%s)\n",ArgIn);

  if (ArgIn.empty()) return 1;
  rangeList_.clear();

  // Set rangeArg
  rangeArg_.assign(ArgIn);
  // Check if ArgIn is a mask expression
  size_t maskcharpos = rangeArg_.find_first_of(":@*");
  if (maskcharpos != std::string::npos) {
    mprinterr("Error: Using a mask expression for range (%s)\n",ArgIn.c_str());
    mprinterr("Error: Ranges should only contain digits, dashes, and commas (e.g. 3-5,8-10)\n");
    return 1;
  } 
  // Split range by comma
  ArgList CommaList(rangeArg_, ",");
  int err = 0;
  while ( err == 0 ) {
    arg = CommaList.GetStringNext();
    if (arg.empty()) break; // Exit the while loop.
    // Then split by dash
    DashList.SetList(arg, "-");
    R[0] = DashList.getNextInteger(-1);
    R[1] = DashList.getNextInteger(-1);
    if (R[0]==-1) {
      mprinterr("Error: Range::SetRange(%s): Range is -1 for %s\n",ArgIn.c_str(), 
                DashList.ArgLine());
      err=1;
      break;
    }
    upper = R[1];
    if (upper==-1) upper=R[0];
    ++upper; // Want up to and including the upper argument
    if ( this->setRange(R[0], upper) )
      mprintf("Warning: Converting %s to range [%i-%i] is not valid.\n",
              ArgIn.c_str(), R[0], R[1]);
  }

  // Do not return an empty list
  if ( err > 0 || rangeList_.empty() ) 
    return 1;
  
  // Sort frames using default comparison
  std::sort(rangeList_.begin(), rangeList_.end());
  //for (it=rangeList_.begin(); it!=rangeList_.end(); it++)
  //  fprintf(stdout,"RangeList= %i\n",*it); 
  // Remove duplicates
  RangeType::const_iterator it = std::unique( rangeList_.begin(), rangeList_.end() );
  rangeList_.resize( it - rangeList_.begin() );

  return 0;
}

// Range::setRange()
/** Given a start and end number, add numbers to range from start to (but
  * not including) end.  
  */
int Range::setRange(int start, int end) {
  //mprintf("DEBUG: setRange called with start=%i end=%i\n", start, end);
  //Check that end is greater than start so that the range list is
  if (end <= start) {
    mprintf("Error: Range::SetRange: end (%i) <= start (%i)\n",end,start);
    return 1;
  }
  for (int range=start; range < end; range++)
    rangeList_.push_back(range);

  return 0;
}

/** Given a start and end number, set up range from start up to but
  * not including end.
  */
int Range::SetRange(int start, int end) {
  rangeList_.clear();
  //mprintf("DEBUG: SetRange called with start=%i end=%i\n", start, end);
  if (setRange(start, end)) return 1;
  // Set range expression
  rangeArg_ = ArrayToRangeExpression( rangeList_, 0 );
  return 0;
}

/** Add number to range. Ensure it is sorted. */
void Range::AddToRange(int num) {
  //mprintf("DEBUG: AddToRange(%i)\n", num);
  //PrintToStdout();
  // Find where num should be
  unsigned int idx = 0;
  while (idx < rangeList_.size()) {
    // Check if num is already present
    if (num == rangeList_[idx]) return;
    if (num < rangeList_[idx]) break;
    idx++;
  }
  //mprintf("DEBUG: %i belongs at index %u\n", num, idx);
  // If it belongs at the end, add it and leave
  if (idx == rangeList_.size()) {
    rangeList_.push_back( num );
    return;
  }
  // Add placeholder to end of the range array
  unsigned int jdx = rangeList_.size();
  rangeList_.push_back(0);
  // Move everything up one
  for (unsigned int i = jdx; i > idx; i--)
    rangeList_[i] = rangeList_[i-1];
  // Put num at idx
  rangeList_[idx] = num;
}

// Range::ShiftBy()
/** Shift all numbers in range by specified value. */
void Range::ShiftBy(int val) {
  for (RangeType::iterator rangeNum = rangeList_.begin(); 
                            rangeNum != rangeList_.end(); ++rangeNum)
    *rangeNum += val;
}

// Range::PrintRange()
/** Print all numbers in the range to a line. Increment by offset. */
void Range::PrintRange(const char* header, int offset) const {
  if (header!=0)
    mprintf("%s",header);
  std::string rangeExp = ArrayToRangeExpression( rangeList_, offset );
  mprintf(" %s", rangeExp.c_str());
/*
  for (std::list<int>::const_iterator it=rangeList_.begin(); it!=rangeList_.end(); it++)
    mprintf(" %i",(*it)+offset);
  //mprintf("\n");
*/
}

/** For debugging, print entire range to stdout. */
void Range::PrintToStdout() const {
  for (RangeType::const_iterator it=rangeList_.begin(); it!=rangeList_.end(); it++)
    mprintf(" %i", *it);
  mprintf("\n");
}

/** \return True if number idx is in the range. */
bool Range::InRange(int idx) const {
  return ( std::find( rangeList_.begin(), rangeList_.end(), idx ) != rangeList_.end() );
}
