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

  mprintf("DEBUG: '%s' %u elts, %u chars wide, %u elts per line, %u lines per block, fsize= %u\n",
          fnameIn.full(), Nelts_, eltWidth_, Ncols_, linesPerBlock_, fsize);

  return BufferedLine::OpenFileRead(fnameIn, fsize);
}
