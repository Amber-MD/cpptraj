#ifndef INC_FILETYPES_H
#define INC_FILETYPES_H
#include "BaseIOtype.h"
#include "ArgList.h"
/** This static class is used to create arrays for allocatable file classes
  * like DataFiles and Trajectories.
  */
class FileTypes {
  public:
    enum OptType { READOPT=0, WRITEOPT };
    typedef int FileFormatType;
    /// Allocator and description for file types.
    /** This array should be constructed with one entry per allocatable
      * type, with the last entry corresponding to UNKNOWN. The UNKNOWN entry
      * should only have a Description, with 0 for all other fields.
      */
    struct AllocToken {
      const char* Description;
      BaseIOtype::HelpType ReadHelp;
      BaseIOtype::HelpType WriteHelp;
      BaseIOtype::AllocatorType Alloc;
    };
    typedef const AllocToken* AllocPtr;
    /// For associating keywords/extensions with file types.
    /** This array can have as many entries as there are recognized keywords/
      * extensions, but must end with the UNKNOWN entry, with the Key and 
      * Extension fields set to 0. When multiple entries exist for a given
      * Allocator, the first one is considered the default.
      */
    struct KeyToken {
      FileFormatType Type;
      const char* Key;
      const char* Extension;
    };
    typedef const KeyToken* KeyPtr;
    /// \return type if file format keyword is present in ArgList.
    static FileFormatType GetFormatFromArg(KeyPtr, ArgList&, FileFormatType);
    /// \return type for given file format keyword.
    static FileFormatType GetFormatFromString(KeyPtr, std::string const&, FileFormatType);
    /// \return default extension for the given type.
    static std::string GetExtensionForType(KeyPtr, FileFormatType);
    /// \return type that matches given extension.
    static FileFormatType GetTypeFromExtension(KeyPtr, std::string const&, FileFormatType);
    /// \return Description for given type.
    static const char* FormatDescription(AllocPtr, FileFormatType);
    /// \return Allocator for given type. MUST BE CAST TO PROPER TYPE.
    static BaseIOtype* AllocIO(AllocPtr, FileFormatType, bool);
    /// List all formats or help for a specific format.
    static void Options(KeyPtr, AllocPtr, FileFormatType, std::string const&, OptType);
  private:
    static std::string FormatKeywords(KeyPtr, FileFormatType);
    static std::string FormatExtensions(KeyPtr, FileFormatType);
    static inline std::string FmtString(KeyPtr, FileFormatType);
};
#endif
