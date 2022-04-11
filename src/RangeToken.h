#ifndef INC_RANGETOKEN_H
#define INC_RANGETOKEN_H
#include <string>
/// Used to split a string of format <num1>[-<num2>]; either number may be negative.
class RangeToken {
  public:
    RangeToken() {}
    /// \return First number string detected by Assign()
    std::string const& Num1() const { return num1_; }
    /// \return Second number string detected by Assign()
    std::string const& Num2() const { return num2_; }
    /// \return Error message if Assign() returns non-zero
    std::string const& Err() const  { return err_; }
    /// \return 0 if range string was converted to numbers, 1 if an error occurred.
    int Assign(std::string const& token) {
      num1_.clear();
      num2_.clear();
      err_.clear();
      bool dash1 = false;
      bool dash2 = false;
      int idx = 0;
      for (std::string::const_iterator it = token.begin(); it != token.end(); ++it)
      {
        if (*it == '-') {
          if (idx == 0) {
            // Initial dash indicating first number negative
            dash1 = true;
            idx = 1;
          } else if (idx == 1) {
            // Should be separating 2 numbers
            if (num1_.empty()) {
              err_.assign("Too many dashes before first number in token '" + token + "'");
              return 1;
            }
            idx = 2;
          } else if ( idx == 2) {
            // Initial dash after separator dash indicating second number negative
            dash2 = true;
            idx = 3;
          } else {
            if (num2_.empty())
              err_.assign("Too many dashes before second number in token '" + token + "'");
            else 
              err_.assign("Trailing dashes in token '" + token + "'");
            return 1;
          }
        } else if (isdigit(*it)) {
          if (idx == 0) {
            idx = 1;
            num1_ += *it;
          } else if (idx == 1) {
            num1_ += *it;
          } else {
            if (idx == 2) idx = 3;
            num2_ += *it;
          }
        }
      }
      if (num1_.empty() && num2_.empty()) {
        err_.assign("No numbers in token '" + token + "'");
        return 1;
      }
      if (dash1) num1_.assign("-" + num1_);
      if (dash2) {
        if (num2_.empty()) {
          err_.assign("Trailing dash after first number in token '" + token + "'");
          return 1;
        }
        num2_.assign("-" + num2_);
      }
      return 0;
    }
  private:
    std::string num1_; ///< First number string
    std::string num2_; ///< Second number string
    std::string err_;  ///< Error message string
};
#endif
