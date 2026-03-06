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
    /// \return true if full path matches
    bool operator==(FileName const& rhs) const { return (fullPathName_ == rhs.fullPathName_); }
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
    /// \return the full file path
    const std::string& Full()      const { return fullPathName_;         }
    /// \return only the base file name
    const std::string& Base()      const { return baseName_;             }
    /// \return only the base file name extension
    const std::string& Ext()       const { return extension_;            }
    /// \return the full file path (const char*)
    const char* full()             const { return fullPathName_.c_str(); }
    /// \return only the base file name (const char*)
    const char* base()             const { return baseName_.c_str();     }
    /// \return only the base file name extension (const char*)
    const char* ext()              const { return extension_.c_str();    }
    /// \return Recognized compression extension if present
    const std::string& Compress()  const { return compressExt_;          }
    /// \return the directory prefix, including the trailing slash
    const std::string& DirPrefix() const { return dirPrefix_;            }
    /// \return the directory prefix, no trailing slash
    std::string DirPrefix_NoSlash() const {
      if (dirPrefix_.empty()) return dirPrefix_;
      return dirPrefix_.substr(0, dirPrefix_.size()-1);
    }
    /// \return true if the file name is not set
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
  /// Search for file names <base>.X given a base name (and debug level) 
  NameArray SearchForReplicas(FileName const&, int);
# ifdef MPI
  /// Base name, trajComm master, ensemble rank, ensemble size, debug
  NameArray SearchForReplicas(FileName const&, bool, int, int, int);
# endif
  /// Print error message corresponding to 'false' value from 'Exists()'
  void ErrorMsg(const char*);
  bool Exists(std::string const&);
  bool Exists(FileName const&);
  int Remove(FileName const&);
  /** Given lowest replica traj filename, split into components for search. */
  class RepName {
  public:
    RepName() : ExtWidth_(0), lowestRepnum_(-1), extChar_('.') {}
    RepName(FileName const&, int);
    bool Error() const { return Prefix_.empty(); }
    /// \return Replica file name for given offset from lowest replica number.
    FileName RepFilename(int) const;
  private:
    std::string Prefix_;      ///< File name up to the numerical extension.
    std::string ReplicaExt_;  ///< Numerical extension.
    std::string CompressExt_; ///< Optional compression extension after numerical extension.
    int ExtWidth_;            ///< Width of the numerical extension. TODO remove
    int lowestRepnum_;        ///< Integer value of numerical extension.
    char extChar_;            ///< Character preceding numerical extension
  };
}
#endif
