#ifndef INC_FILENAME_H
#define INC_FILENAME_H
#include <string>
/// Class to hold file name, extension, etc.
class FileName {
  public:
    FileName() {}
    /// CONSTRUCTOR - Take previously set FileName and a tag.
    FileName(FileName const&, std::string const&);
    FileName(const FileName&);
    FileName& operator=(const FileName&);
    /// Set file name and extensions; only name is known.
    int SetFileName(std::string const&);
    /// Set file name and extension; perform tilde expansion.
    int SetFileNameWithExpansion(std::string const&);
    /// Set file name and extensions; name and compressed status are known.
    int SetFileName(std::string const&, bool);
    /// Set file tag.
    void SetTag(std::string const& tIn) { tag_ = tIn; }
    /// Clear FileName
    void clear();
    /// \return true if input tag/name is a match for this FileName.
    bool IsMatch(std::string const&) const;

    const std::string& Full()     const { return fullPathName_;         }
    const std::string& Base()     const { return baseName_;             }
    const char* full()            const { return fullPathName_.c_str(); }
    const char* base()            const { return baseName_.c_str();     }
    const std::string& Tag()      const { return tag_;                  }
    const std::string& Ext()      const { return extension_;            }
    const std::string& Compress() const { return compressExt_;          }
    bool empty()                  const { return fullPathName_.empty(); }
  private:
    enum CompressStatus { UNKNOWN=0, YES, NO };
    std::string fullPathName_;
    std::string baseName_;
    std::string extension_;
    std::string compressExt_;
    std::string tag_;

    int SetFileName(std::string const&, CompressStatus);
};
#endif
