#ifndef INC_FILELIST_H
#define INC_FILELIST_H
#include <vector>
#include <string>
class FileList {
  public:
    FileList();
    virtual ~FileList();

    inline void SetDebug(int dIn) {
      debug_ = dIn;
    }
    int FindTag(std::string&);
    int FindName(std::string&);
  protected:
    int debug_;
  private:
    std::vector<std::string> tags_;
    std::vector<std::string> fnames_;
    std::vector<std::string> basenames_;
};
#endif
