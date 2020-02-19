#ifndef INC_FILETYPE_H
#define INC_FILETYPE_H
#include <string>
namespace Cpptraj {
namespace File {
  enum Type { T_FILE = 0, T_DIR, T_UNKNOWN };
  /// \return File type
  Type IdType(std::string const&);
}
}
#endif
