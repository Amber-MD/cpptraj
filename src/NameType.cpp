#include <algorithm> // std::fill, std::copy
#include "NameType.h"

NameType::NameType()
{
  std::fill(c_array_, c_array_ + NameSize_-1, ' ');
  c_array_[NameSize_-1] = '\0';
}

NameType::NameType(const NameType &rhs)
{
  std::copy(rhs.c_array_, rhs.c_array_ + NameSize_, c_array_);
}

NameType::NameType(const char *rhs)
{
  const char *ptr = rhs;
  for (unsigned int j = 0; j < NameSize_; j++) {
    c_array_[j] = *ptr;
    if (*ptr=='\0') break;
    ++ptr;
  }
  FormatName();
}

NameType::NameType(std::string const& str)
{
  unsigned int ns1 = NameSize_ - 1;
  unsigned int strend = (unsigned int)str.size();
  if (strend > ns1)
    strend = ns1;
  for (unsigned int j = 0; j < strend; j++) 
    c_array_[j] = str[j];
  c_array_[strend] = '\0';
  FormatName();
}
 
NameType &NameType::operator=(const NameType &rhs) {
  if (&rhs==this) return *this;
  std::copy(rhs.c_array_, rhs.c_array_ + NameSize_, c_array_);
  return *this;
}

// NameType::ToBuffer()
/** For interfacing with old C stuff. Only set 1st 4 chars. */
void NameType::ToBuffer(char *buffer) const {
  buffer[0] = c_array_[0];
  buffer[1] = c_array_[1];
  buffer[2] = c_array_[2];
  buffer[3] = c_array_[3];
  buffer[4] = '\0';
}

bool NameType::Match(NameType const& maskName) const { 
  int c = 0;
  for (unsigned int m = 0; m < NameSize_-1; m++) {
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

bool NameType::operator==(const NameType &rhs) const {
  for (unsigned int idx = 0; idx < NameSize_; idx++) {
    if (c_array_[idx] != rhs.c_array_[idx]) return false;
    // If we are here, chars at idx (including if null) must be equal.
    if (c_array_[idx] == '\0') break;
  }
  return true;
}

bool NameType::operator==(const char *rhs) const {
  NameType tmp(rhs);
  return (*this == tmp);
}

bool NameType::operator!=(const NameType &rhs) const {
  for (unsigned int idx = 0; idx < NameSize_; idx++) {
    if (c_array_[idx] != rhs.c_array_[idx]) return true;
    // If we are here, chars at idx (including if null) must be equal.
    if (c_array_[idx] == '\0') break;
  }
  return false;
}

bool NameType::operator!=(const char *rhs) const {
  NameType tmp(rhs);
  return (*this != tmp);
}

char NameType::operator[](int idx) const {
  if (idx < 0 || idx >= (int)NameSize_) return '\0';
  return c_array_[idx];
}

std::string NameType::Truncated() const {
  unsigned int i = 0;
  for (; i != NameSize_; i++)
    if (c_array_[i] == ' ' || c_array_[i] == '\0') break;
  return std::string( c_array_, c_array_+i );
}

/** \return Non-space length of name. */
int NameType::len() const {
  unsigned int i = 0;
  for (; i != NameSize_; i++)
    if (c_array_[i] == ' ' || c_array_[i] == '\0')
      return (int)i;
  return (int)i;
}

/** Replace asterisks with a single quote */
void NameType::ReplaceAsterisk() {
  for (unsigned int idx = 0; idx < NameSize_; idx++)
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
}
