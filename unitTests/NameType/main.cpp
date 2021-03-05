// Unit test for NameType class
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>
#include "NameType.h"

static const int Err(const char* msg) {
  fprintf(stderr, "Error: %s\n", msg);
  return 1;
}

struct reverseSort {
  bool operator() (const NameType& lhs, const NameType& rhs) const {
    return ( lhs > rhs );
  }
} reverseSortObj;

int main() {
  // This unit test is primarily for the matching.
  NameType name1("CT11");
  // Copy construction
  NameType tmp("HA");
  NameType name2(tmp);
  // Test assignment
  NameType name3;
  name3 = name2;
  // Test conversion from strings
  std::string str4("HA2");
  NameType name4( str4 );
  std::string str5("ZN");
  NameType name5( str5 );

  // Copy construction
  std::vector<NameType> Nvector;
  Nvector.push_back( name3 );
  Nvector.push_back( name5 );
  Nvector.push_back( name4 );
  Nvector.push_back( name1 );
  Nvector.push_back( name2 );

  // Test that assignment worked.
  bool n1n2_eq = (name2 == name3);
  if (!n1n2_eq) return Err("'==' operator failed after assignment.\n");
  bool n1n2_ne = (name2 != name3);
  if (n1n2_ne) return Err("'!=' operator failed after assignment.\n");

  n1n2_eq = (name2 == "HA");
  if (!n1n2_eq) return Err("'==' const char* operator failed after assignment.\n");
  n1n2_ne = (name2 != "HA");
  if (n1n2_ne) return Err("'!=' const char* operator failed after assignment.\n");

  // Sort. Tests the < operator.
  std::sort( Nvector.begin(), Nvector.end() );

  // Test that copy construction/sort worked.
  if ( (name1 != Nvector[0]) ||
       (name2 != Nvector[1]) ||
       (name3 != Nvector[2]) ||
       (name4 != Nvector[3]) ||
       (name5 != Nvector[4]) )
    return Err("Sort/copy construction failed (order is wrong).\n");
 
  // Reverse sort. Tests the > operator.
  std::sort( Nvector.begin(), Nvector.end(), reverseSortObj );

  // Test that reverse sort worked.
  if ( (name5 != Nvector[0]) ||
       (name4 != Nvector[1]) ||
       (name3 != Nvector[2]) ||
       (name2 != Nvector[3]) ||
       (name1 != Nvector[4]) )
    return Err("Reverse sort failed (order is wrong).\n");

  // Write to stdout
  for (std::vector<NameType>::const_iterator it = Nvector.begin(); it != Nvector.end(); ++it)
    printf(" '%s'", *(*it));
  printf("\n");

  // Test that an out of range index is properly caught
  char testchar1 = name1[-1];
  char testchar2 = name1[999999];
  if (testchar1 != '\0' || testchar2 != '\0')
    return Err("[] operator did not return null char for out of range index.\n");

  // Test length function 
  if (name1.len() != 4 ||
      name2.len() != 2 ||
      name3.len() != 2 ||
      name4.len() != 3 ||
      name5.len() != 2)
    return Err("len() function failed.\n");

  // Test single wildcard matching for name1 'CT11'
  bool swc1 = name1.Match("?T11");
  bool swc2 = name1.Match("C?11");
  bool swc3 = name1.Match("CT?1");
  bool swc4 = name1.Match("CT1?");
  bool swc5 = name1.Match("?T1?");
  if (!swc1 ||
      !swc2 ||
      !swc3 ||
      !swc4 ||
      !swc5)
    return Err("Single wildcard matching failed.\n");

  // Test wildcard matching for name1 'CT11'
  bool wc1 = name1.Match("C*");
  bool wc2 = name1.Match("C*1");
  bool wc3 = name1.Match("CT*");
  if (!wc1 ||
      !wc2 ||
      !wc3)
    return Err("Wildcard matching failed.\n");

  // Test mixed wildcard matching
  bool mwc1 = name1.Match("*1?");
  if (!mwc1)
    return Err("Mixed wildcard matching failed.\n");

  // Check that large name produces a warning
  NameType large1("ThisIsALargeName");
  NameType large2( std::string("ThisIsALargeString") );

  // This should not produce a warning
  NameType char5("12345");

  // Test matching with asterisk in name.
  NameType name6("O5*1");
  bool a1wc = name6.Match("O5*");    // Match
  bool a2wc = name6.Match("O5\\*1"); // Match
  bool a3wc = name6.Match("O5\\*?"); // Match
  if (!a1wc ||
      !a2wc ||
      !a3wc)
    return Err("Matching name with asterisk.\n");
  NameType name7("O541");
  bool a4wc = name7.Match("O5*");    // Match
  bool a5wc = name7.Match("O5\\*1"); // No match
  bool a6wc = name7.Match("O5\\*?"); // No Match
  if (!a4wc ||
       a5wc ||
       a6wc)
    return Err("Matching name without asterisk.\n");

  return 0;
}
