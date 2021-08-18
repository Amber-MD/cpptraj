#ifndef INC_TEXTBLOCKBUFFER_H
#define INC_TEXTBLOCKBUFFER_H
#include "BufferedLine.h"
/// Used to read in regular blocks of text from file
class TextBlockBuffer : private BufferedLine {
  public:
    /// CONSTRUCTOR
    TextBlockBuffer();
    /// Open file for reading and set up block buffer for # elts, elt width, elts per line, additional bytes
    int OpenFileRead(FileName const&, unsigned int, unsigned int, unsigned int, unsigned int);
    /// Set up text block for given # elts, elt width, elts per line
    int SetupTextBlock(unsigned int, unsigned int, unsigned int);
    /// Read block into double array
    int BlockToDoubles(double*);
    // Members of BufferedLine that should be public
    using BufferedLine::Filename;
    using BufferedLine::CloseFile;
    using BufferedLine::Line;
    using BufferedLine::CurrentLine;
    using BufferedLine::Nremaining;
    using BufferedLine::SetDebug;
  private:
    using BufferedLine::IsDos;
    using BufferedLine::Debug;

    unsigned int Nelts_;         ///< Number of elements in a block
    unsigned int Ncols_;         ///< Max number of elements on a line
    unsigned int eltWidth_;      ///< Width of each element
    unsigned int linesPerBlock_; ///< Number of lines in a block
};
#endif
