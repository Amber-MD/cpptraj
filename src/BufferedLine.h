#ifndef INC_BUFFEREDLINE_H
#define INC_BUFFEREDLINE_H
#include <vector>
#include "CpptrajFile.h"
/// Used to buffer text files that will be read line-by-line
class BufferedLine : private CpptrajFile {
  public:
    BufferedLine();
    ~BufferedLine();

    const char* Line();
    int TokenizeLine(const char*);
    const char* NextToken();
    /// Open file for reading, set up buffer.
    int OpenFileRead( std::string const& fname ) {
      if ( OpenRead( fname ) ) return 1;
      return ResetBuffer();
    }
    int LineNumber()     const { return nline_;  }
    const char* Buffer() const { return buffer_; }
    // Members of CpptrajFile that should be public
    using CpptrajFile::Filename;
    using CpptrajFile::CloseFile;
    using CpptrajFile::GetLine;
  private:
    int ResetBuffer();
    static const size_t DEFAULT_BUFFERSIZE = 16384;

    char* buffer_;         ///< Character buffer
    char* bufferPosition_; ///< Position in buffer/start of current line.
    /// Array of pointers to beginning and ends of tokens in current line. 
    std::vector<char*> tokens_;
    size_t tokenidx_;      ///< Current position in tokens array
    char saveChar_;        ///< Saved last char of current token
    char* lineEnd_;        ///< End of current line in buffer
    char endChar_;         ///< Character that was at *lineend
    char* endBuffer_;      ///< End position of buffer
    size_t nline_;         ///< Current line number.
};
#endif
