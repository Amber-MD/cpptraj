#include <algorithm> // std::min
#include "TextBlockBuffer.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
TextBlockBuffer::TextBlockBuffer() :
  Nelts_(0),
  Ncols_(0),
  eltWidth_(0),
  linesPerBlock_(0)
{}

/** Open file and set up for reading regular text block. */
int TextBlockBuffer::OpenFileRead(FileName const& fnameIn, unsigned int nelts,
                                  unsigned int eltwidth, unsigned int eltsPerLine,
                                  unsigned int additionalBytes)
{
  Nelts_ = nelts;
  Ncols_ = eltsPerLine;
  eltWidth_ = eltwidth;

  linesPerBlock_ = Nelts_ / Ncols_;
  if ((Nelts_ % Ncols_) > 0)
    ++linesPerBlock_;

  // Determine size of each frame;
  unsigned int nEndChars;
  if (IsDos())
    nEndChars = 2;
  else
    nEndChars = 1;

  unsigned int fsize = (Nelts_ * eltWidth_) + (linesPerBlock_ * nEndChars) + additionalBytes;
  if (Debug() > 0)
    mprintf("DEBUG: '%s' %u elts, %u chars wide, %u elts per line, %u lines per block, fsize= %u\n",
            fnameIn.full(), Nelts_, eltWidth_, Ncols_, linesPerBlock_, fsize);

  return BufferedLine::OpenFileRead(fnameIn, fsize);
}

/** Set up size of text block. */
int TextBlockBuffer::SetupTextBlock(unsigned int nelts,
                                    unsigned int eltwidth, unsigned int eltsPerLine)
{
  Nelts_ = nelts;
  Ncols_ = eltsPerLine;
  eltWidth_ = eltwidth;

  linesPerBlock_ = Nelts_ / Ncols_;
  if ((Nelts_ % Ncols_) > 0)
    ++linesPerBlock_;

  if (Debug() > 0)
    mprintf("DEBUG: '%s' %u elts, %u chars wide, %u elts per line, %u lines per block\n",
            Filename().full(), Nelts_, eltWidth_, Ncols_, linesPerBlock_);

  return 0;
}

/** Read block into double array.
  * \return Number of elements actually read.
  */
int TextBlockBuffer::BlockToDoubles(double* darray) {
  unsigned int elt = 0;
  // Loop over lines in block
  char* ptrend = 0;
  for (unsigned int line = 0; line != linesPerBlock_; line++)
  {
    Line();
    char* ptr = BufferPosition();
    if (*ptr == '\0' || ptr == 0) return (int)elt;
    unsigned int nRemaining = Nelts_ - elt;
    unsigned int maxcol = std::min(Ncols_, nRemaining);
    // Loop over columns in line
    for (unsigned int col = 0; col < maxcol; col++)
    {
      ptrend = ptr + eltWidth_;
      char lastchar = *ptrend;
      *ptrend = '\0';
      //mprintf("DEBUG: %12u : '%s'\n", elt+1, ptr);
      darray[elt++] = atof(ptr);
      *ptrend = lastchar;
      ptr = ptrend;
    }
  }
  // Update buffer position with final position
  if (ptrend != 0) SetBufferPosition(ptrend);
  return (int)elt;
}
