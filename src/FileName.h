#ifndef INC_FILENAME_H
#define INC_FILENAME_H
#include <string>
#include <vector>
/// Class to hold file name, extension, etc. TODO add to File namespace?
class FileName {
  public:
    FileName() {}
    FileName(std::string const& s) { SetFileName(s); }
    FileName(const char* s) { SetFileName( std::string(s) ); }
    FileName(const FileName&);
    FileName& operator=(const FileName&);
    /// Set file name and extensions, perform expansion as necessary.
    int SetFileName(std::string const&);
    /// Set file name, no expansions.
    int SetFileName_NoExpansion(std::string const&);
    /// Append given string to file name but do not change extension info.
    int Append(std::string const&); //TODO this can probably replace AppendNumber some places
    /// \return Filename with given string appended - do not change extension info.
    FileName AppendFileName(std::string const&) const;
    /// \return Filename with given string prepended to base file name.
    FileName PrependFileName(std::string const&) const;
    /// \return Filename with given string prepended to extension
    FileName PrependExt(std::string const&) const;
    /// Clear FileName
    void clear();
    /// \return true if string matches full or base file name.
    bool MatchFullOrBase(std::string const&) const;

    const std::string& Full()      const { return fullPathName_;         }
    const std::string& Base()      const { return baseName_;             }
    const std::string& Ext()       const { return extension_;            }
    const char* full()             const { return fullPathName_.c_str(); }
    const char* base()             const { return baseName_.c_str();     }
    const char* ext()              const { return extension_.c_str();    }
    const std::string& Compress()  const { return compressExt_;          }
    const std::string& DirPrefix() const { return dirPrefix_;            }
    bool empty()                   const { return fullPathName_.empty(); }
  private:
    std::string fullPathName_;
    std::string baseName_;
    std::string extension_;
    std::string compressExt_;
    std::string dirPrefix_;
};
/// This namespace contains useful file-related routines.
namespace File {
  typedef std::vector<FileName> NameArray;
  /// Expand given expression to array of file names
  NameArray ExpandToFilenames(std::string const&);
  /// Print error message corresponding to 'false' value from 'Exists()'
  void ErrorMsg(const char*);
  bool Exists(std::string const&);
  bool Exists(FileName const&);
}
#endif
