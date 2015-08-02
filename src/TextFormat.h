#ifndef INC_TEXTFORMAT_H
#define INC_TEXTFORMAT_H
#include <cstddef> // size_t
#include <string>
/// Hold printf-like text format string.
class TextFormat {
  public:
    /// Formats: f, E, g, i, s
    enum FmtType { DOUBLE = 0, SCIENTIFIC, GDOUBLE, INTEGER, STRING };
    /// CONSTRUCTOR - default 8.3f format
    TextFormat() : type_(DOUBLE), width_(8), precision_(3), nelements_(1),
                   leftAlign_(false), isLong_(false) {}
    /// Set format string of type given width and precision, optionally append
    void SetFormatString(FmtType, int, int, bool);
    /// Set format string of type given width and precision
    void SetFormatString(FmtType t, int w, int p) {
      SetFormatString(t, w, p, false);
    }
    /// Append double format string for size, min, step, default width and precision 
    void SetCoordFormat(size_t, double, double, int, int);
    /// \return Format string.
    const char* fmt() const { return fmt_.c_str(); }
  private:
    static char TypeChar_[]; ///< Hold printf format chars for each type.

    std::string fmt_; ///< Hold format string.
    FmtType type_;    ///< Format type.
    int width_;       ///< Width of each element.
    int precision_;   ///< precision of each element.
    int nelements_;   ///< Number of elements.
    bool leftAlign_;  ///< Alignment of each element.
    bool isLong_;     ///< If true format is long, needed for double reads only.
};
#endif
