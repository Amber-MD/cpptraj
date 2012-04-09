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

    void AddNames(char *, const char*, std::string &);
    void AddNames(char *, std::string &, std::string &);
    void AddFilename(char *);
    int FindName(char*); // TODO: Make obsolete
    int FindName(std::string);
    // TODO: Place these in definition
    std::string &Tag(int idx) {
      return tags_[idx];
    }
    std::string &Name(int idx) {
      return basenames_[idx];
    }
    bool HasNames() {
      if (!tags_.empty() || !fnames_.empty() || !basenames_.empty())
        return true;
      return false;
    }
  protected:
    int debug_;
  private:
    std::vector<std::string> tags_;
    std::vector<std::string> fnames_;
    std::vector<std::string> basenames_;
};
#endif
