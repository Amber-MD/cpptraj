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
    mprintf("Warning: Name truncation detected: Name='%s' vs Original='%s'\n", c_array_, rhs);
  }
}

/** Initialize NameType with given buffer. */
NameType::NameType(const char *rhs)
{
  if (rhs != 0)
    Assign( rhs );
  else
    c_array_[0] = '\0';
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
    if (maskName.c_array_[m] == '\0' && c_array_[c] == '\0') {
      // At end of names, all must have matched up until now.
      break;
    } else if (maskName.c_array_[m] == '\0' || c_array_[c] == '\0') {
      // One name has ended without the other ending; the only way this
      // can match is if one of the names ends with *. Otherwise mismatch.
      if (maskName.c_array_[m] == '*' || c_array_[c] == '*')
        return true;
      else
        return false;
    }
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

/** NOTE: NameType is always "Truncated" now (no extra whitespace),
  *       but the function name is retained to avoid a lot of
  *       rewriting for the sake of rewriting.
  */
std::string NameType::Truncated() const {
  return std::string( c_array_ );
}

/** \return name formatted; left-aligned with given minimal width. If
  *         necessary spaces are added to make minimum width.
  */
std::string NameType::Formatted(int minWidth) const {
  std::string nameOut;
  nameOut.reserve(ArraySize_);
  nameOut.assign( c_array_ );
  int nblank = minWidth - (int)nameOut.size();
  if (nblank > 0)
    nameOut.append(nblank, ' ');
  return nameOut;
}

/** \return Non-space length of name. */
int NameType::len() const {
  const char* ptr = c_array_;
  while (*ptr != '\0') ++ptr;
  return (int)(ptr - c_array_);
}
