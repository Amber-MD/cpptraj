#include <cstring> // strchr 
#include "BufferedLine.h"
#include "CpptrajStdio.h"

BufferedLine::BufferedLine() :
  currentBufSize_(DEFAULT_BUFFERSIZE),
  buffer_(0),
  bufferPosition_(0),
  tokenidx_(0),
  saveChar_(0),
  lineEnd_(0),
  endBuffer_(0),
  nline_(0)
{}

BufferedLine::~BufferedLine() {
  if (buffer_!=0) delete[] buffer_;
}

// BufferedLine::ResetBuffer()
int BufferedLine::ResetBuffer() {
  if (buffer_!=0) delete[] buffer_;
  buffer_ = new char[ currentBufSize_ ];
  bufferPosition_ = buffer_;              // Point to beginning of buffer.
  endBuffer_ = buffer_ + currentBufSize_; // Point to 1 past end of buffer.
  lineEnd_ = endBuffer_;                  // Indicates buffer needs to be filled.
  nline_ = 0;
  return 0;
} 

/** Open file for reading with given buffer size. */
int BufferedLine::OpenFileRead(FileName const& fnameIn, size_t bufsize) {
  currentBufSize_ = bufsize;
  if (OpenRead( fnameIn )) return 1;
  return ResetBuffer();
}

// BufferedLine::Line()
const char* BufferedLine::Line() {
  bufferPosition_ = lineEnd_;
  //mprintf("DEBUG: currentBufSize= %zu  buffer= %zu  bufferPosition= %zu"
  //        "  lineEnd= %zu  endBuffer= %zu\n",
  //        currentBufSize_, buffer_, bufferPosition_ - buffer_, 
  //        lineEnd_ - buffer_, endBuffer_ - buffer_);
  // Search for the end of the next line in buffer_
  while (lineEnd_ <= endBuffer_) {
    // Fill buffer if needed
    if (lineEnd_ == endBuffer_) {
      // Shift this line to the beginning and fill remainder of buffer.
      size_t bufferRemainder = endBuffer_ - bufferPosition_;
      //mprintf("DEBUG: bufferRemainder %zu\n", bufferRemainder);
      if (bufferRemainder == currentBufSize_) break;
      std::copy(bufferPosition_, bufferPosition_ + bufferRemainder, buffer_);
      int Nread = Read(buffer_ + bufferRemainder, currentBufSize_ - bufferRemainder);
      //mprintf("DEBUG: Attempted read of %zu bytes, actually read %i bytes (remainder %zu).\n",
      //        currentBufSize_ - bufferRemainder, Nread, bufferRemainder);
      // Set buffer position
      bufferPosition_ = buffer_;
      lineEnd_ = buffer_ + bufferRemainder;
      endBuffer_ = lineEnd_ + (size_t)Nread;
      if (Nread < 1) {
        // If nothing was read and no previous buffer, set null.
        if (bufferRemainder == 0) buffer_[0] = '\0';
        return 0;
      }
    }
    if ( *lineEnd_ == '\n') { // End of the line. Replace with null char.
      *(lineEnd_++) = '\0';
      ++nline_;
      return bufferPosition_;
    }
    ++lineEnd_;
  }

  // If we are here there was no newline. Increase the buffer size and try again.
  bool reallocate = true;
  while (reallocate) {
    size_t new_buf_size = currentBufSize_ * 2;
    char* new_buffer = new char[ new_buf_size ];
    endBuffer_ = new_buffer + new_buf_size;
    std::copy( buffer_, buffer_ + currentBufSize_, new_buffer );
    bufferPosition_ = new_buffer + (bufferPosition_ - buffer_);
    delete[] buffer_;
    buffer_ = new_buffer;
    //mprintf("DEBUG: Buffer has been reallocated: new size %zu\n", new_buf_size);
    // Try to fill the remainder of the new buffer.
    int Nread = Read(buffer_ + currentBufSize_, new_buf_size - currentBufSize_);
    //mprintf("DEBUG: Attempted additional read of %zu bytes, actually read %i bytes.\n",
    //        new_buf_size - currentBufSize_, Nread);
    if (Nread < 1) {
      // No additonal read possible. Warn and return the current buffer.
      mprintf("Warning: No newline in file.\n");
      lineEnd_ = buffer_ + currentBufSize_;
      *lineEnd_ = '\0';
      currentBufSize_ = new_buf_size;
      ++nline_;
      return bufferPosition_;
    }
    endBuffer_ = buffer_ + currentBufSize_ + (size_t)Nread; 
    // Search for newline in additonal read
    lineEnd_ = buffer_ + currentBufSize_;
    //mprintf("DEBUG: endBuffer %zu  lineEnd %zu\n", endBuffer_ - buffer_, lineEnd_ - buffer_);
    currentBufSize_ = new_buf_size;
    while (lineEnd_ != endBuffer_) {
      if (*lineEnd_ == '\n') {
        *(lineEnd_++) = '\0';
        ++nline_;
        return bufferPosition_;
      }
      ++lineEnd_;
    }
  }
  return 0; 
}

/** Separate the current line into tokens delimited by given chars. 
  * \return Number of tokens.
  */
int BufferedLine::TokenizeLine(const char* separator) {
  if ( separator == 0 ) return 0;
  char* linechar = bufferPosition_;
  bool inToken = false;
  tokens_.clear();
  // NOTE: Just check for newline?
  while ( *linechar != '\0' && *linechar != '\n' ) {
    if (!inToken) { // Not in token.
      if ( strchr( separator, *linechar ) == 0 ) {
        tokens_.push_back( linechar); // Pointer to beginning of token.
        inToken = true;
      }
    } else {        // In a token.
      if ( strchr( separator, *linechar ) != 0 ) {
        tokens_.push_back( linechar ); // Pointer to end of token.
        inToken = false;
      }
    }
    ++linechar;
  }
  // If inToken is still true point to linechar as the last token
  if (inToken)
    tokens_.push_back(linechar);
  tokenidx_ = 0; 
  /*mprintf("DBG: Tokenize: Line=[%s]\n", bufferPosition_);
  mprintf("\t%i Tokens:\n", ntokens);
  for (unsigned int t = 0; t < ntokens; ++t)
    mprintf("\t\t%u %c\n",t, *tokens_[t]);*/
  return (int)(tokens_.size() / 2);
}

/** Following a call to TokenizeLine return char* to next token in the
  * line. Token is generated by inserting a null char at the current
  * token end position; this char is restored once NextToken is called
  * again.
  */
const char* BufferedLine::NextToken() {
  if (tokenidx_ == tokens_.size()) return 0;
  char* tokenptr = tokens_[tokenidx_];
  if (tokenidx_ != 0)
    *(tokens_[tokenidx_-1]) = saveChar_;
  char* nextptr = tokens_[tokenidx_+1];
  saveChar_ = *nextptr;
  *nextptr = '\0';
  tokenidx_ += 2;
  return tokenptr;
}
