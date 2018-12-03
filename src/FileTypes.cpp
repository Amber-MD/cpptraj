#include <set>
#include <algorithm> // std::max
#include <cstring> // strlen
#include "FileTypes.h"
#include "CpptrajStdio.h"
// FileTypes::GetFormatFromArg()
FileTypes::FileFormatType FileTypes::GetFormatFromArg(KeyPtr begin, ArgList& argIn,
                                                      FileFormatType def)
{
  for (KeyPtr token = begin; token->Key != 0; ++token)
    if (argIn.hasKey( token->Key )) return token->Type;
  return def;
}
// FileTypes::GetFormatFromString()
FileTypes::FileFormatType FileTypes::GetFormatFromString(KeyPtr begin, std::string const& fmt,
                                                         FileFormatType def)
{
  for (KeyPtr token = begin; token->Key != 0; ++token)
    if ( fmt.compare( token->Key )==0 ) return token->Type;
  return def;
}
// FileTypes::GetExtensionForType()
std::string FileTypes::GetExtensionForType(KeyPtr begin, FileFormatType typeIn) {
  for (KeyPtr token = begin; token->Extension != 0; ++token)
    if ( token->Type == typeIn )
      return std::string( token->Extension );
  return std::string();
}
// FileTypes::GetTypeFromExtension()
FileTypes::FileFormatType FileTypes::GetTypeFromExtension(KeyPtr begin, std::string const& extIn,
                                                          FileFormatType def)
{
  for (KeyPtr token = begin; token->Extension != 0; ++token)
    if ( extIn.compare( token->Extension ) == 0 ) return token->Type;
  return def;
}
// -----------------------------------------------------------------------------
// FileTypes::FormatDescription()
const char* FileTypes::FormatDescription(AllocPtr allocArray, FileFormatType typeIn) {
  return allocArray[ typeIn ].Description;
}
// FileTypes::AllocIO()
/** \param silent If true do not print error message if allocator does not exist.
  */
BaseIOtype* FileTypes::AllocIO(AllocPtr allocArray, FileFormatType typeIn, bool silent) {
  if (allocArray[typeIn].Alloc == 0) {
    if (!silent)
      mprinterr("Error: CPPTRAJ was compiled without support for %s files.\n",
                allocArray[typeIn].Description);
    return 0;
  }
  return allocArray[typeIn].Alloc();
}
// FileTypes::FormatKeywords()
std::string FileTypes::FormatKeywords(KeyPtr begin, FileFormatType ftype) {
  std::string keywords;
  std::set<std::string> Keys;
  for (KeyPtr token = begin; token->Key != 0; ++token)
    if (token->Type == ftype) Keys.insert( std::string(token->Key) );
  if (!Keys.empty()) {
    keywords.assign("Keywords:");
    for (std::set<std::string>::const_iterator key = Keys.begin();
                                               key != Keys.end(); ++key)
    {
      if (key != Keys.begin()) keywords.append(",");
      keywords.append(" " + *key);
    }
  }
  return keywords;
}
// FileTypes::FormatExtensions()
std::string FileTypes::FormatExtensions(KeyPtr begin, FileFormatType ftype) {
  std::string extensions;
  std::set<std::string> Exts;
  for (KeyPtr token = begin; token->Extension != 0; ++token)
    if (token->Type == ftype) Exts.insert( std::string(token->Extension) );
  if (!Exts.empty()) {
    extensions.assign("Extensions:");
    for (std::set<std::string>::const_iterator ext = Exts.begin();
                                               ext != Exts.end(); ++ext)
    {
      if (ext != Exts.begin()) extensions.append(",");
      extensions.append(" '" + *ext + "'");
    }
  }
  return extensions;
}

// FileTypes::Options()
/** \return 1 if all options listed, 0 if specific option was listed. */
void FileTypes::Options(KeyPtr begin, AllocPtr allocArray, FileFormatType UNK,
                            std::string const& fkey, OptType otype)
{
  if (fkey.empty()) {
    // Everything
    // For nice formatting, get size of largest string
    unsigned int maxsize = 0;
    for (int i = 0; i < UNK; i++)
      maxsize = std::max(maxsize, (unsigned int)strlen(allocArray[i].Description));
    for (int i = 0; i < UNK; i++) {
      mprintf("      %*s:", maxsize, allocArray[i].Description);
      std::string fmtKeywords   = FormatKeywords(begin, i);
      std::string fmtExtensions = FormatExtensions(begin, i);
      if (!fmtExtensions.empty()) fmtKeywords.append(";");
      mprintf(" %s %s\n", fmtKeywords.c_str(), fmtExtensions.c_str());
    }
  } else {
    // Specific format
    FileFormatType ft = GetFormatFromString(begin, fkey, UNK);
    if (ft == UNK)
      mprintf("    Invalid format specifier: %s\n", fkey.c_str());
    else {
      int i = (int)ft;
      std::string fmtKeywords   = FormatKeywords(begin, i);
      std::string fmtExtensions = FormatExtensions(begin, i);
      if (!fmtExtensions.empty()) fmtKeywords.append(",");
      mprintf("    Options for %s: %s %s\n", allocArray[i].Description,
              fmtKeywords.c_str(), fmtExtensions.c_str());
      switch (otype) {
        case READOPT:
          if (allocArray[i].ReadHelp != 0) allocArray[i].ReadHelp(); break;
        case WRITEOPT:
          if (allocArray[i].WriteHelp != 0) allocArray[i].WriteHelp(); break;
      }
    }
  }
}

// FileTypes::WriteOptions()
void FileTypes::WriteOptions(KeyPtr begin, AllocPtr allocArray, FileFormatType UNK) {
  for (int i = 0; i < UNK; i++) {
    std::string fmtExtensions = FormatExtensions(begin, i);
    std::string fmtKeywords =  FormatKeywords(begin, i);
    if (allocArray[i].WriteHelp || !fmtExtensions.empty() || !fmtKeywords.empty()) {
      mprintf("    Options for %s:", allocArray[i].Description);
      if (!fmtKeywords.empty()) mprintf(" %s,", fmtKeywords.c_str());
      if (!fmtExtensions.empty()) mprintf(" %s", fmtExtensions.c_str()); 
      mprintf("\n");
      if (allocArray[i].WriteHelp != 0) allocArray[i].WriteHelp();
    }
  }
}
