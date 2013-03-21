#ifndef INC_FILELIST_H
#define INC_FILELIST_H
#include <vector>
#include <string>
class FileList {
  public:
    FileList();
    virtual ~FileList();
    void Clear();
    inline void SetDebug(int dIn)          { debug_ = dIn;           }
    std::string const& Tag(int idx)  const { return tags_[idx];      }
    std::string const& Name(int idx) const { return basenames_[idx]; }

    void AddNameWithTag(std::string const&, std::string const&, std::string const&);
    void AddFilename(std::string const&);
    int FindName(std::string const&) const;
    bool HasNames() const {return (!tags_.empty() || !fnames_.empty() || !basenames_.empty());}
  protected:
    int debug_;
  private:
    std::vector<std::string> tags_;
    std::vector<std::string> fnames_;
    std::vector<std::string> basenames_;
};
#endif
