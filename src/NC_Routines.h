#ifndef INC_NC_ROUTINES_H
#define INC_NC_ROUTINES_H
#ifdef BINTRAJ
#include <string>
namespace NC {
  /// \return true if given code is error and print message, false otherwise.
  bool CheckErr(int);
  /// \return Text for attribute of given variable ID.
  std::string GetAttrText(int, int, const char*);
  /// \return Text for given global attribute.
  std::string GetAttrText(int, const char*);
  /// \return dimension ID of given attribute and set dimension length.
  int GetDimInfo(int, const char*, unsigned int&);
  // FIXME This version here for backwards compatibility.
  int GetDimInfo(int, const char*, int&);
  /// Print debug info to STDOUT
  void Debug(int);
}
#endif /* BINTRAJ */
#endif
