// Name.cpp
// Collection of routines for manipulating NAME type vars
#include "Name.h"
#include <cctype>

// =============================================================================
// NOTE: The following NAME_ routines may eventually replace the NAME
//       datatype. Not doing this now since ptraj routines are still
//       being used and it is easier and faster to just send them
//       a pointer to POD type.
const int NAME_::NAME_SIZE = 6;

NAME_::NAME_() {
  c_array[0]=' ';
  c_array[1]=' ';
  c_array[2]=' ';
  c_array[3]=' ';
  c_array[4]=' ';
  c_array[5]='\0';
}

NAME_::NAME_(const NAME_ &rhs) {
  c_array[0] = rhs.c_array[0];
  c_array[1] = rhs.c_array[1];
  c_array[2] = rhs.c_array[2];
  c_array[3] = rhs.c_array[3];
  c_array[4] = rhs.c_array[4];
  // 5 is always NULL
}

NAME_ &NAME_::operator=(const NAME_ &rhs) {
  if (&rhs==this) return *this;
  c_array[0] = rhs.c_array[0];
  c_array[1] = rhs.c_array[1];
  c_array[2] = rhs.c_array[2];
  c_array[3] = rhs.c_array[3];
  c_array[4] = rhs.c_array[4];
  // 5 is always NULL
  return *this;
}

char *NAME_::Assign(char *bufferIn, int width) {
  int max;
  char *ptr = bufferIn;
  if (width > NAME_SIZE)
    max = NAME_SIZE;
  else
    max = width;
  for (int j = 0; j < max; j++) {
    c_array[j] = *ptr;
    ++ptr;
  }
  c_array[max-1]='\0';
  return ptr;
}

const char *NAME_::operator*() const {
  return c_array;
}
// =============================================================================

// PadWithSpaces()
/** For consistency with Amber names, replace any NULL in the first 4 chars
  * with spaces.
  */
void PadWithSpaces(NAME name) {
  if (name[0]=='\0') {
    name[0]=' ';
    name[1]=' ';
    name[2]=' ';
    name[3]=' ';
  } else if (name[1]=='\0') {
    name[1]=' ';
    name[2]=' ';
    name[3]=' ';
  } else if (name[2]=='\0') {
    name[2]=' ';
    name[3]=' ';
  } else if (name[3]=='\0') {
    name[3]=' ';
  }
  name[4]='\0';
}

// TrimName()
/** Remove leading whitespace from NAME */
void TrimName(NAME NameIn) {
  if        (NameIn[0]!=' ') { // No leading whitespace
    return; 
  } else if (NameIn[1]!=' ') { // [_XXX]
    NameIn[0]=NameIn[1];
    NameIn[1]=NameIn[2];
    NameIn[2]=NameIn[3];
    NameIn[3]=' ';
  } else if (NameIn[2]!=' ') { // [__XX]
    NameIn[0]=NameIn[2];
    NameIn[1]=NameIn[3];
    NameIn[2]=' ';
    NameIn[3]=' ';
  } else if (NameIn[3]!=' ') { // [___X]
    NameIn[0]=NameIn[3];
    NameIn[1]=' ';
    NameIn[2]=' ';
    NameIn[3]=' ';
  } 
  // Otherwise Blank, no trim needed
  return;
} 

// WrapName()
/** Move leading characters that are digits to the back of the string until 
  * first char is alphabetic.
  */
void WrapName(NAME NameIn) {
  int i;
  int blank=-1;
  int numalpha=0;
  char digit;
  //fprintf(stderr,"NameIn = [%s]\n",NameIn);
  // First trim leading whitespace if necessary
  // NOTE: Should already be done
  //TrimName(NameIn);
  // Check for last blank char and alpha chars
  for (i=0; i<4; i++) {
    if (!isdigit(NameIn[i])) {
      // If the first char is already alpha exit
      if (i==0) return;
      numalpha++;
    }
    if (blank<0 && NameIn[i]==' ') blank=i;
  }
  if (blank==-1) blank=4;
  // No alpha chars, return
  if (numalpha==0) return;
  // No chars (____), just return
  if (blank==0) return;
  // 1 char (X___), also return 
  if (blank==1) return;
  // Otherwise 2 char (XX__), 3 char (XXX_), or 4 char (XXXX)
  while (isdigit(NameIn[0])) {
    digit=NameIn[0];
    for (i=1; i<blank; i++) NameIn[i-1]=NameIn[i];
    NameIn[blank-1]=digit;
  }
  //fprintf(stderr,"NameOut= [%s]\n",NameIn);
}

// ReplaceAsterisk()
/** Given a string of length 5 (4chars + 1 NULL) change any asterisk (*) to
  * prime ('). In cpptraj asterisks are considered reserved characters for
  * atom masks.
  */
void ReplaceAsterisk(NAME NameIn) {
  if (NameIn[0]=='*') NameIn[0]='\'';
  if (NameIn[1]=='*') NameIn[1]='\'';
  if (NameIn[2]=='*') NameIn[2]='\'';
  if (NameIn[3]=='*') NameIn[3]='\'';
  return;
}

