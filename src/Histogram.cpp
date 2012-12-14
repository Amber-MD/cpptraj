#include <cmath> // log
#include <stdexcept>
#include "Histogram.h"
#include "CpptrajStdio.h"
#include "CpptrajFile.h"
// Histogram

// CONSTRUCTOR
Histogram::Histogram() : 
  debug_(0),
  isCircular_(0),
  hasCycled_(false)
{
  // DataSet-specific vars
  width_ = 12;
  precision_ = 4;
  dType_ = HIST;
  SetDataSetFormat(false);
}

// Histogram::SetDebug()
void Histogram::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_>0) mprintf("Histogram DEBUG LEVEL SET TO %i\n", debug_);
}

Dimension& Histogram::operator[](int idx) {
  if (idx < 0 || idx >= (int)dimensions_.size()) {
    mprinterr("Error: Histogram: Dimension %i out of bounds.\n",idx);
    throw std::out_of_range("Histogram[]");
  }
  return dimensions_[idx];
}

// Histogram::AddDimension()
/** Add a dimension to the histogram with the given min, max, step,
  * and number of bins.
  * NOTE: NO ERROR CHECKING IS PERFORMED.
  */
int Histogram::AddDimension(Dimension &dim) {

  dimensions_.push_back( dim );

  if (debug_>0) {
    mprintf("ADDED DIM:");
    dimensions_.back().PrintDim();
  }

  // NOTE: Should the following be its own routine?
  // Recalculate offsets for all dimensions starting at farthest coord. This 
  // follows column major ordering. 
  int offset = 1; 
  for (std::vector<Dimension>::reverse_iterator rd = dimensions_.rbegin();
                                                rd != dimensions_.rend(); ++rd)
  {
    if (debug_>0) mprintf("\tHistogram: %s offset is %i\n",(*rd).c_str(), offset);
    (*rd).SetOffset( offset );
    offset *= (*rd).Bins(); 
  }
  // offset should now be equal to the total number of bins across all dimensions
  if (debug_>0) mprintf("\tHistogram: Total Bins = %i\n",offset);
  // Allocate space for bins and set to 0
  Bins_.resize(offset, 0);

  // Update number of DataSet dimensions
  dim_ = (int)dimensions_.size();

  return 0;
}

// Histogram::BinData()
/** Given an array of data with the same dimensionality as the histogram, bin
  * the data.
  */
int Histogram::BinData(std::vector<double>& DataIn) {
  int index=0;
  std::vector<double>::iterator data = DataIn.begin();

  // DEBUG
  //mprintf("Binning [");
  //for (int i=0; i<numDimension; i++) mprintf("%lf,",Data[i]);
  //mprintf("]\n");
  // Loop over defined dimensions. 
  // Calculate an index into Bins based on precalcd offsets for dimensions.
  // Populate bin.
  if (debug_>1) mprintf("\t{");
  for (std::vector<Dimension>::iterator dim = dimensions_.begin();
                                        dim != dimensions_.end(); ++dim)
  {
    // Check if Data is out of bounds for this dimension 
    if (*data > (*dim).Max() || *data < (*dim).Min()) {
      index = -1;
      break;
    }
    // Calculate index (idx) for this particular dimension 
    double coord = *data - (*dim).Min();
    coord = coord / (*dim).Step();
    int idx = (int)coord;
    if (debug_>1) mprintf(" [%s:%lf (%i)],",(*dim).c_str(), *data, idx);

    // Calculate overall index in Bins, offset has already been calc.
    index += (idx * (*dim).Offset());
    ++data;
  }

  // If index was successfully calculated, populate bin 
  if (index>-1 && index < (int)Bins_.size()) {
    if (debug_>1) mprintf(" |index=%i",index);
    ++Bins_[index];
  } else {
    mprintf("\tWarning: Coordinates out of bounds %i {",index);
    for (data = DataIn.begin(); data != DataIn.end(); ++data) 
      mprintf(" %f", *data);
    mprintf(" }\n");
  }

  if (debug_>1) mprintf("}\n");
  return 0;
}

// -----------------------------------------------------------------------------
// Histogram::BinStart()
/** Set current bin to 0 and initialize indices. If isCircularIn is true the
  * bin indices will wrap in each dimension.
  */
void Histogram::BinStart(bool circularIn) {
  BinIndices_.clear();
  if (circularIn) {
    isCircular_ = 1;
    BinIndices_.resize( dimensions_.size(), -1 );
  } else {
    isCircular_ = 0;
    BinIndices_.resize( dimensions_.size(),  0 );
  }
  //mprintf("BinStart: Dimension %zu, circular is %i\n",BinIndices_.size(), isCircular_);
  //mprintf("D1=%i  D2=%i\n",dimensions_[0].Bins(), dimensions_[1].Bins());
}

// Histogram::PrintBins()
/** This routine will print out the contents of the Bin array in column-major 
  * order. A counter is used to calculate the appropriate coordinate indices 
  * for the array. The count array works like this:
  *   Dim0index Dim1index ... DimNindex
  * DimNindex is considered the highest order dimension and is always 
  * incremented. When DimNindex reaches the number of bins for that Dim
  * it cycles, and the next highest order dim is checked, down to the lowest
  * order dimension.
  * If circular specified, wrapping will occur so that data from indices -1 
  * and N+1 (corresponding to N and 0) are printed out as well.
  * --- CURRENTLY NOT IMPLEMENTED ---
  * If binType is 0, the Bins array will be printed. If binType is 1, the
  * landscape will be printed.
  * if SD is not null, the standard deviation array will also be printed.
  */
void Histogram::PrintBins(const char *outfilename, bool circularIn, bool gnuplot) {
  CpptrajFile outfile;
  if (outfile.SetupWrite(outfilename, debug_)) return;
  if (outfile.OpenFile()) return;

  mprintf("\tHistogram: Writing standard histogram file %s\n",outfilename);

  BinStart( circularIn );

  if (gnuplot) {
    if (dimensions_.size() == 2) 
      outfile.Printf("set pm3d map\nsplot \"-\" with pm3d title \"%s\"\n",outfilename);
    else if (dimensions_.size() == 1)
      outfile.Printf("plot \"-\"\n",outfilename);
  }

  if (debug_>0) {
    if (circularIn)
      mprintf("\t\tPrinting %zu bins in circular fashion.\n", Bins_.size());
    else
      mprintf("\t\tPrinting %zu bins.\n",Bins_.size());
  }

  bool loop=true;
  while (loop) {
    int index = BinIndicesToIndex();
    // If we dont care about zero bins or bin pop > 0, output
    for (unsigned int i=0; i < dimensions_.size(); ++i)
      outfile.Printf("%lf ", 
                      ((double)BinIndices_[i]*dimensions_[i].Step()) + dimensions_[i].Min() );
    outfile.Printf("%lf\n",Bins_[index]);

    loop = IncrementBinIndices( );
    // If gnuplot, print extra space when highest order coord cycles
    if (gnuplot && hasCycled_) 
        outfile.Printf("\n");
  }
  if (gnuplot && dimensions_.size() < 3)
    outfile.Printf("end\npause -1\n"); 
  outfile.CloseFile();
}

// Histogram::BinIndicesToIndex()
/** Calculate index based on BinIndices, converting wrapped indices to 
  * actual indices.
  */
int Histogram::BinIndicesToIndex() {
  int idx;
  int index = 0;
  std::vector<int>::iterator count = BinIndices_.begin();
  for (std::vector<Dimension>::iterator dim = dimensions_.begin(); 
                                        dim != dimensions_.end(); ++dim)
  {
    //mprinterr(" %i",count[n]);
    if (*count == -1) 
      idx = (*dim).Bins() - 1;
    else if (*count == (*dim).Bins()) 
      idx = 0;
    else 
      idx = *count;
    index += (idx * (*dim).Offset());
    ++count;
  }
  //mprinterr(" = %i\n",index);
  return index;
}

// Histogram::IncrementBinIndices()
/** \return true if there are more bins to process.
  * \return false if there are no more bins.
  */
bool Histogram::IncrementBinIndices() {
  //mprintf("DEBUG0\t\t\tCoord0=%i Coord1=%i\n",BinIndices_[0],BinIndices_[1]);
  // Increment highest order coord.
  std::vector<int>::reverse_iterator rcount = BinIndices_.rbegin(); 
  ++(*rcount);
  std::vector<Dimension>::reverse_iterator rdim = dimensions_.rbegin();
  // Check if highest order coord has Cycled.
  if (*rcount == (*rdim).Bins() + isCircular_)
    hasCycled_ = true;
  else
    hasCycled_ = false;
  // Increment other coords if necessary
  for (; rdim != dimensions_.rend() - 1; ++rdim) {
    if (*rcount == (*rdim).Bins()+isCircular_) {
      (*rcount) = 0 - isCircular_; 
      ++rcount;
      ++(*rcount);
    }
  }
  //mprintf("DEBUG1\t\t\tCoord0=%i Coord1=%i\n",BinIndices_[0],BinIndices_[1]);
  // If the lowest order coord is at lowest dimensions #bins, done.
  if (BinIndices_[0] == dimensions_[0].Bins() + isCircular_) 
    return false;

  return true;
}

void Histogram::WriteBuffer(CpptrajFile& cbuffer, int frame) {
  if (frame < 0 || frame >= (int)Bins_.size())
    cbuffer.Printf(data_format_, 0.0);
  else
    cbuffer.Printf(data_format_, Bins_[frame]);
}

void Histogram::Write2D(CpptrajFile& outfile, int x, int y) {
  if ( x < 0 || x >= dimensions_[0].Bins())
    outfile.Printf(data_format_, 0.0);
  else if ( y < 0 || y >= dimensions_[1].Bins() )
    outfile.Printf(data_format_, 0.0);
  else {
    // TODO: Use BinIndicesToIndex
    int index = (x * dimensions_[0].Offset()) + (y * dimensions_[1].Offset());
    outfile.Printf(data_format_, Bins_[index]);
  }
}

void Histogram::GetDimensions(std::vector<int>& vIn) {
  vIn.clear();
  vIn.reserve( dimensions_.size() );
  for (std::vector<Dimension>::iterator dim = dimensions_.begin(); 
                                        dim != dimensions_.end(); ++dim)
  {
    vIn.push_back( (*dim).Bins() );
  }
}

// -----------------------------------------------------------------------------
// Histogram::CalcFreeE()
/** Calculate free energy based on bin populations.  */
int Histogram::CalcFreeE(double T, int refbin) {
  mprintf("\tHistogram: Calculating free E at %lf K.\n",T);
  // TODO: Make Kb a constant
  double KT = (-.001985886596 * T);

  // Find most populated bin for G=0
  std::vector<double>::iterator bin = Bins_.begin();
  double binmax = *bin;
  ++bin;
  for (; bin != Bins_.end(); ++bin) 
    if (*bin > binmax) 
      binmax = *bin;
  mprintf("\t           Bins max is %.0lf\n",binmax);
  if (binmax==0) {
    mprinterr("Histogram: Cannot calc free E, no bins populated!\n");
    return 1;
  }

  // If requested, set up reference bin other than max
  if (refbin>-1) {
    if (Bins_[refbin] > 0) {
      binmax = Bins_[refbin];
      mprintf("\t           Calculating free E w.r.t bin %i, population %lf\n",refbin,binmax);
    } else
      mprintf("Warning: Reference bin %i has no population. Using %lf\n",refbin,binmax);
  }

  // Set artificial ceiling for bins with 0 population. Make it equivalent
  // to a bin population of 0.5. 
  double temp = 0.5 / binmax;
  double ceiling = log(temp) * KT;
  //ceiling*=1.10;
  mprintf("\t           Artificial ceiling is %lf kcal/mol.\n",ceiling);

  // Calculate free E based on populations
  for (bin = Bins_.begin(); bin != Bins_.end(); ++bin) {
    temp = *bin;               // Store Bin population in temp
    if (temp>0) {
      temp /= binmax;          // bin/max
      *bin = log(temp) * KT;   // -R*T*ln(bin/max)
    } else
      *bin = ceiling;          // Artificial ceiling for 0 pop bins.
  }

  return 0;
}

// Histogram::Normalize()
/** Normalize bins so that sum over all bins is 1.0 */
int Histogram::Normalize() {
  double sum = 0.0;
  mprintf("\tHistogram: Normalizing bin populations to 1.0\n");
  for (std::vector<double>::iterator bin = Bins_.begin(); bin != Bins_.end(); ++bin)
    sum += *bin;
  mprintf("\t           Sum over all bins is %lf\n",sum);
  if (sum == 0.0) {
    mprinterr("Error: Histogram::Normalize: Sum over bin populations is 0.0\n");
    return 1;
  }
  for (std::vector<double>::iterator bin = Bins_.begin(); bin != Bins_.end(); ++bin)
    *bin /= sum;
  return 0;
}

