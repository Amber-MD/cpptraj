#ifndef INC_TEXTFORMAT_H
#define INC_TEXTFORMAT_H
#include <cstddef> // size_t
#include <string>
/// Hold printf-like text format string.
class TextFormat {
  public:
    /// Formats: f, E, g, i, s
    enum FmtType { DOUBLE = 0, SCIENTIFIC, GDOUBLE, INTEGER, STRING };
    /// CONSTRUCTOR - default 8.3f format - no setup
    TextFormat() : type_(DOUBLE), width_(8), precision_(3), nelements_(1),
                   leftAlign_(false), isLong_(false) {}
    /// CONSTRUCTOR - type, width - no setup
    TextFormat(FmtType t, int w) : type_(t), width_(w), precision_(0), nelements_(1),
               leftAlign_(false), isLong_(false) {}
    /// CONSTRUCTOR - type, width, leftalign
    TextFormat(FmtType t, int w, bool l) : type_(t), width_(w), precision_(0), nelements_(1),
               leftAlign_(l), isLong_(false)
      { SetFormatString(type_, width_, precision_); }
    /// CONSTRUCTOR - For output coords, size, min, step, width, precision
    TextFormat(size_t z, double m, double s, int w, int p) : type_(DOUBLE),
               width_(w), precision_(p), nelements_(1), leftAlign_(false), isLong_(false)
      { SetCoordFormat( z, m, s, w, p ); }
    /// Set format string of type with given width and precision
    void SetFormatString(FmtType, int, int); // TODO do not take any args
    /// Set format string with current type, width, and precision, optional left align
    void SetFormatString(bool l) {
      leftAlign_ = l;
      SetFormatString(type_, width_, precision_);
    }
    /// Set double format string for size, min, step, default width and precision 
    void SetCoordFormat(size_t, double, double, int, int);
    /// Set width
    void SetWidth(int w) { width_ = w; }
    /// Set precision
    void SetPrecision(int p) { precision_ = p; }
    /// \return pointer to format string.
    const char* fmt() const { return fmt_.c_str(); }
    int Width()       const { return width_;       }
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
