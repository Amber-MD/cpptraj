#ifndef INC_FILELIST_H
#define INC_FILELIST_H
#include <vector>
#include "FileName.h"
class FileList {
  public:
    FileList();
    virtual ~FileList();
    void Clear();
    inline void SetDebug(int dIn)          { debug_ = dIn;        }
    std::string const& Tag(int idx)  const { return tags_[idx];   }
    FileName const& Name(int idx)    const { return fnames_[idx]; }

    void AddNameWithTag(FileName const&, std::string const&);
    void AddFilename(FileName const&);
    int FindName(std::string const&) const;
    bool HasNames() const {return (!tags_.empty() || !fnames_.empty());}
  protected:
    int debug_;
  private:
    std::vector<std::string> tags_;
    std::vector<FileName> fnames_;
};
#endif
