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

    void AddNames(char *, std::string &, std::string &);
    int FindName(char*); // TODO: Make obsolete
    int FindName(std::string&);
  protected:
    int debug_;
  private:
    std::vector<std::string> tags_;
    std::vector<std::string> fnames_;
    std::vector<std::string> basenames_;
};
#endif
