#ifndef INC_TEXTFORMAT_H
#define INC_TEXTFORMAT_H
#include <cstddef> // size_t
#include <string>
/// Hold printf-like text format string.
class TextFormat {
  public:
    /// Formats: f, E, g, i, s
    enum FmtType { DOUBLE = 0, SCIENTIFIC, GDOUBLE, INTEGER, STRING };
    /// Alignment type: Right, Left, Right with leading space.
    enum AlignType { RIGHT = 0, LEFT, LEADING_SPACE };
    /// CONSTRUCTOR - default 8.3f format
    TextFormat() : type_(DOUBLE), width_(8), precision_(3), nelements_(1),
                   colwidth_(0), align_(RIGHT), isLong_(false) { SetFormatString(); }
    /// CONSTRUCTOR - Type only
    TextFormat(FmtType t) : type_(t), width_(0), precision_(-1), nelements_(1),
               colwidth_(0), align_(RIGHT), isLong_(false) { SetFormatString(); }
    // TODO trim constructors
    /// CONSTRUCTOR - type, width
    TextFormat(FmtType t, int w) : type_(t), width_(w), precision_(0), nelements_(1),
               colwidth_(0), align_(RIGHT), isLong_(false) { SetFormatString(); }
    /// CONSTRUCTOR - type, width, alignment
    TextFormat(FmtType t, int w, AlignType a) : type_(t), width_(w), precision_(0), nelements_(1),
               colwidth_(0), align_(a), isLong_(false) { SetFormatString(); }
    /// CONSTRUCTOR - type, width, precision
    TextFormat(FmtType t, int w, int p) : type_(t), width_(w), precision_(p), nelements_(1),
               colwidth_(0), align_(RIGHT), isLong_(false) { SetFormatString(); }
    /// CONSTRUCTOR - type, width, precision, nelements
    TextFormat(FmtType t, int w, int p, int n) : type_(t), width_(w), precision_(p), nelements_(n),
               colwidth_(0), align_(RIGHT), isLong_(false) { SetFormatString(); }
    /// CONSTRUCTOR - For output coords, size, min, step, width, precision
    TextFormat(size_t z, double m, double s, int w, int p) : type_(DOUBLE), width_(w),
               precision_(p), nelements_(1), colwidth_(0), align_(RIGHT), isLong_(false)
      { SetCoordFormat( z, m, s, w, p ); }
    /// Set double format string for size, min, step, default width and precision 
    void SetCoordFormat(size_t, double, double, int, int);
    /// Set format string with new alignment.
    void SetFormatAlign(AlignType t) { align_ = t; SetFormatString(); }
    /// Set format string with new width
    void SetFormatWidth(int w) { width_ = w; SetFormatString(); }
    /// Set format string with new width and precision
    void SetFormatWidthPrecision(int w, int p) { width_ = w; precision_ = p; SetFormatString(); }
    /// Set width only - do not reset format string (for DataSet_string).
    void SetWidth(int w) { width_ = w; }
    /// \return pointer to format string.
    const char* fmt()        const { return fmt_.c_str(); }
    std::string const& Fmt() const { return fmt_;         }
    int Width()              const { return width_;       }
    int Precision()          const { return precision_;   }
    int ColumnWidth()        const { return colwidth_;    }
  private:
    /// Set format string for current type, width, precision, etc 
    void SetFormatString();

    static char TypeChar_[]; ///< Hold printf format chars for each type.

    std::string fmt_;  ///< Hold format string.
    FmtType type_;     ///< Format type.
    int width_;        ///< Width of each element.
    int precision_;    ///< precision of each element.
    int nelements_;    ///< Number of elements.
    int colwidth_;     ///< Total width required for format string.
    AlignType align_;  ///< Alignment of each element.
    bool isLong_;      ///< If true format is long, needed for double reads only.
};
#endif
