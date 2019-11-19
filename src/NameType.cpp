#include <cctype>    // isspace
#include <algorithm> // std::copy
#include "NameType.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
NameType::NameType()
{
  c_array_[0] = '\0';
}

/// COPY CONSTRUCTOR
NameType::NameType(const NameType &rhs)
{
  std::copy(rhs.c_array_, rhs.c_array_ + ArraySize_, c_array_);
}

/** Assign incoming buffer (up to ArraySize_) to this NameType. Ignore whitespace. */
void NameType::Assign( const char* rhs ) {
  const char* ptr = rhs;
  unsigned int j = 0;
  while (j < ArraySize_) {
    if (*ptr == '\0') break;
    if (!isspace(*ptr))
      c_array_[j++] = *ptr;
    ptr++;
  }
  // Detect truncation of input.
  if (j < ArraySize_)
    // No truncation
    c_array_[j] = '\0';
  else {
    c_array_[ArraySize_-1] = '\0';
    mprintf("Warning: Name truncation detected: Name='%s' vs Original'%s'\n", c_array_, rhs);
  }
}

/** Initialize NameType with given buffer. */
NameType::NameType(const char *rhs)
{
  if (rhs != 0)
    Assign( rhs );
  else
    c_array_[0] = '\0';
  FormatName();
}

/** Initialize NameType with given string. */
NameType::NameType(std::string const& str)
{
  if (!str.empty())
    Assign(str.c_str());
  else
    c_array_[0] = '\0';
/*
  unsigned int ns1 = ArraySize_ - 1;
  unsigned int strend = (unsigned int)str.size();
  if (strend > ns1)
    strend = ns1;
  for (unsigned int j = 0; j < strend; j++) 
    c_array_[j] = str[j];
  c_array_[strend] = '\0';*/
  FormatName();
}

/// ASSIGNMENT 
NameType &NameType::operator=(const NameType &rhs) {
  if (&rhs==this) return *this;
  std::copy(rhs.c_array_, rhs.c_array_ + ArraySize_, c_array_);
  return *this;
}

/** Copy NameType to buffer. For interfacing with old C stuff.*/
void NameType::ToBuffer(char *buffer) const {
  if (buffer == 0) return;
  unsigned int idx = 0;
  for (const char* ptr = c_array_; *ptr != '\0'; ++ptr)
    buffer[idx++] = *ptr;
  buffer[idx] = '\0';
}

/** See if this NameType matches the given NameType.
  * \param maskName Name to match; may include single '?' or multiple '*' char wildcard(s)
  */
bool NameType::Match(NameType const& maskName) const { 
  int c = 0;
  for (unsigned int m = 0; m < ArraySize_-1; m++) {
    if (maskName.c_array_[m] == '\0' && c_array_[c] == ' ')
      // At end of mask and whitespace in name: OK
      break;
    if (maskName.c_array_[m] == '\\') { 
      // Backslash: match literal next char in mask
      ++m;
      if (maskName.c_array_[m] != c_array_[c])
        return false;
    } else if (maskName.c_array_[m] == '*') { 
      // Mask wildcard: instant match
      return true;
    } else if (maskName.c_array_[m] != '?' && 
               maskName.c_array_[m] != c_array_[c]) { 
      // Not mask single wildcard and mismatch
      return false;
    }
    //mprintf("(%c,%c)",maskName.c_array_[m],c_array_[c]);
    ++c;
  }
  return true;
}

/** \return True only if incoming NameType is an exact match. No wildcards. */
bool NameType::operator==(const NameType &rhs) const {
  for (unsigned int idx = 0; idx < ArraySize_; idx++) {
    if (c_array_[idx] != rhs.c_array_[idx]) return false;
    // If we are here, chars at idx (including if null) must be equal.
    if (c_array_[idx] == '\0') break;
  }
  return true;
}

/** \return True only if incoming NameType is an exact match. No wildcards. */
bool NameType::operator==(const char *rhs) const {
  NameType tmp(rhs);
  return (*this == tmp);
}

/** \return True only if incoming NameType does not match. No wildcards. */
bool NameType::operator!=(const NameType &rhs) const {
  for (unsigned int idx = 0; idx < ArraySize_; idx++) {
    if (c_array_[idx] != rhs.c_array_[idx]) return true;
    // If we are here, chars at idx (including if null) must be equal.
    if (c_array_[idx] == '\0') break;
  }
  return false;
}

/** \return True only if incoming NameType does not match. No wildcards. */
bool NameType::operator!=(const char *rhs) const {
  NameType tmp(rhs);
  return (*this != tmp);
}

/** \return Character at given position, or null if position is out of range. */
char NameType::operator[](int idx) const {
  if (idx < 0 || idx >= (int)ArraySize_) return '\0';
  return c_array_[idx];
}

std::string NameType::Truncated() const {
  unsigned int i = 0;
  for (; i != ArraySize_; i++)
    if (c_array_[i] == ' ' || c_array_[i] == '\0') break;
  return std::string( c_array_, c_array_+i );
}

/** \return Non-space length of name. */
int NameType::len() const {
  unsigned int i = 0;
  for (; i != ArraySize_; i++)
    if (c_array_[i] == ' ' || c_array_[i] == '\0')
      return (int)i;
  return (int)i;
}

/** Replace asterisks with a single quote */
void NameType::ReplaceAsterisk() {
  for (unsigned int idx = 0; idx < ArraySize_; idx++)
  {
    if (c_array_[idx] == '\0') break;
    if (c_array_[idx] == '*') c_array_[idx]='\'';
  }
}

// NameType::FormatName()
/** For consistency with Amber names, replace any null in the first 4 chars
  * with spaces. Remove any leading whitespace.
  */
void NameType::FormatName() 
{
  // Remove leading whitespace.
  // Find index of first non-whitespace (blank) char.
  unsigned int nonWSidx = 0;
  while (c_array_[nonWSidx] == ' ')
    ++nonWSidx;
  if (nonWSidx > 0) {
    unsigned int idx = 0;
    for (; nonWSidx < ArraySize_; ++nonWSidx, ++idx) {
      c_array_[idx] = c_array_[nonWSidx];
      if (c_array_[idx] == '\0') break;
    }
  }
  // Ensure at least 4 chars long.
  if (c_array_[0]=='\0') { // 0 chars
    c_array_[0]=' ';
    c_array_[1]=' ';
    c_array_[2]=' ';
    c_array_[3]=' ';
    c_array_[4]='\0';
  } else if (c_array_[1]=='\0') { // 1 char
    c_array_[1]=' ';
    c_array_[2]=' ';
    c_array_[3]=' ';
    c_array_[4]='\0';
  } else if (c_array_[2]=='\0') { // 2 chars
    c_array_[2]=' ';
    c_array_[3]=' ';
    c_array_[4]='\0';
  } else if (c_array_[3]=='\0') { // 3 chars
    c_array_[3]=' ';
    c_array_[4]='\0';
  }
/*
  // Remove leading whitespace.
  if (c_array_[0]==' ') { // Some leading whitespace
    if (c_array_[1]!=' ') {        // [_XXX]
      c_array_[0]=c_array_[1];
      c_array_[1]=c_array_[2];
      c_array_[2]=c_array_[3];
      c_array_[3]=' ';
    } else if (c_array_[2]!=' ') { // [__XX]
      c_array_[0]=c_array_[2];
      c_array_[1]=c_array_[3];
      c_array_[2]=' ';
      c_array_[3]=' ';
    } else if (c_array_[3]!=' ') { // [___X]
      c_array_[0]=c_array_[3];
      c_array_[1]=' ';
      c_array_[2]=' ';
      c_array_[3]=' ';
    }
  }
*/
}
