#ifndef INC_REMOTE_H
#define INC_REMOTE_H
#include <string>
namespace Cpptraj {
/// Class for downloading remote files.
class Remote {
  public:
    /// CONSTRUCTOR
    Remote();
    /// CONSTRUCTOR - Take remote base URL
    Remote(std::string const&);
    /// Set whether to overwrite files or not.
    void SetOverwrite(bool);
    /// Set debug level
    void SetDebug(int);
    /// Download file to output file assuming URL is remote directory
    int DownloadFile(std::string const&, std::string const&) const;
    /// Download file assuming URL is remote directory
    int DownloadFile(std::string const&) const;
  private:
    /// Set remote download command
    int setRemoteDownloadCommand();

    bool overwrite_;           ///< If true overwrite any existing files
    int debug_;
    std::string url_;          ///< Remote base URL
    static std::string cmd_;   ///< Command to call to download remote files
    static std::string oflag_; ///< Command output file flag
};
}
#endif
