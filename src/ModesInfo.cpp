#include <cstdio> //sscanf
#include <cstring> //strncmp
#include "ModesInfo.h"
#include "CpptrajStdio.h"
#include "CpptrajFile.h"

const size_t ModesInfo::BUFSIZE_ = 1024;

// CONSTRUCTOR
ModesInfo::ModesInfo() :
  type_(MT_UNKNOWN),
  source_(MS_UNKNOWN),
  navgelem_(0),
  avg_(0),
  nvect_(0),
  nvectelem_(0),
  freq_(0),
  evec_(0)
{
  dType_ = MODES;
}

// CONSTRUCTOR
ModesInfo::ModesInfo(modesType tIn, modesSource sIn, std::string& nameIn) :
  type_(tIn),
  source_(sIn),
  navgelem_(0),
  avg_(0),
  nvect_(0),
  nvectelem_(0),
  freq_(0),
  evec_(0)
{
  // DataSet
  name_ = nameIn;
  dType_ = MODES;
}

int ModesInfo::SetNavgElem(int mask1tot) {
  // TODO: This is already calcd in MatrixType Nelt.
  if (type_ == MT_DIST || type_ == MT_IDEA || type_ == MT_IRED)
    navgelem_ = mask1tot;
  else if (type_ == MT_DISTCOVAR)
    navgelem_ = mask1tot * (mask1tot - 1) / 2;
  else // CORREL, COVAR, MWCOVAR
    navgelem_ = 3 * mask1tot;
  //mprintf("DEBUG: ModesInfo::navgelem = %i\n",navgelem_);
  return navgelem_;
}

// DESTRUCTOR
ModesInfo::~ModesInfo() {
  if (avg_!=0) delete[] avg_;
  if (freq_!=0) delete[] freq_;
  if (evec_!=0) delete[] evec_;
}

// ModesInfo::ReadEvecFile()
int ModesInfo::ReadEvecFile(std::string& modesfile, int ibeg, int iend) {
  CpptrajFile infile;
  char buffer[BUFSIZE_];

  if (infile.SetupRead( modesfile.c_str(), 0 )) return 1;
  if (infile.OpenFile()) return 1;
  // Read title line
  if (infile.IO->Gets(buffer, BUFSIZE_)!=0) {
    mprinterr("Error: ReadEvecFile(): error while reading title (%s)\n",infile.Name());
    return 1;
  }
  source_ = MS_FILE;

  std::string title(buffer);
  if (title.find("DISTCOVAR")!=std::string::npos)
    type_ = MT_DISTCOVAR;
  else if (title.find("MWCOVAR")!=std::string::npos)
    type_ = MT_MWCOVAR;
  else if (title.find("COVAR")!=std::string::npos)
    type_ = MT_COVAR;
  else if (title.find("DIST")!=std::string::npos)
    type_ = MT_DIST;
  else if (title.find("CORREL")!=std::string::npos)
    type_ = MT_CORREL;
  else if (title.find("IDEA")!=std::string::npos)
    type_ = MT_IDEA;
  else if (title.find("IRED")!=std::string::npos)
    type_ = MT_IRED;
  else {
    // For compatibility with quasih and nmode output
    mprintf("Warning: ReadEvecFile(): Unrecognized type [%s]\n", title.c_str());
    mprintf("         Assuming MWCOVAR.\n");
    type_ = MT_MWCOVAR;
  }

  // Read number of coords for avg and evec
  if ( infile.IO->Gets(buffer, BUFSIZE_)!=0) {
    mprinterr("Error: ReadEvecFile(): error while reading number of atoms (%s)\n",infile.Name());
    return 1;
  }
  int nvals = sscanf(buffer, "%i %i", &navgelem_, &nvectelem_);
  if (nvals == 0) {
    mprinterr("Error: ReadEvecFile(): sscanf on coords failed (%s)\n",infile.Name());
    return 1;
  } else if (nvals == 1) {
    mprintf("Warning: ReadEvecFile(): No value for nvectelem found in %s, assuming it is navgelem\n",
            infile.Name());
    nvectelem_ = navgelem_;
  }

  // Allocate memory for avg, freq, evec
  if (avg_!=0) delete[] avg_;
  if (freq_!=0) delete[] freq_;
  if (evec_!=0) delete[] evec_;
  avg_ = new double[ navgelem_ ];
  // TODO: Check bounds
  freq_ = new double[ iend - ibeg + 1 ];
  evec_ = new double[ (iend - ibeg + 1) * nvectelem_ ];

  // Read avg coordinates
  int ncoords = navgelem_;
  int nlines = ncoords / 7;
  if (ncoords > 0 && (ncoords % 7) != 0)
    ++nlines;
  int nent = 0;
  double tmpval[7];
  for (int i = 0; i < nlines; ++i) {
    if (infile.IO->Gets(buffer, BUFSIZE_)!=0) {
      mprinterr("Error: ReadEvecFile(): error while reading avg coords (%s)\n",infile.Name());
      return 1;
    }
    nvals = sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf", tmpval, tmpval+1, tmpval+2,
                   tmpval+3, tmpval+4, tmpval+5, tmpval+6);
    for (int j = 0; j < nvals; j++)
      avg_[nent++] = tmpval[j];
  }

  // Read eigenvectors
  ncoords = nvectelem_;
  nlines = ncoords / 7;
  if (ncoords > 0 && (ncoords % 7) != 0)
    ++nlines;
  nvect_ = 0;
  int nno = 0;
  while ( infile.IO->Gets(buffer, BUFSIZE_)!=0 ) { // This should read in ' ****'
    if (strncmp(buffer," ****", 5)!=0) {
      mprinterr("Error: ReadEvecFile(): When reading eigenvectors, expected ' ****',\n");
      mprinterr("       got %s [%s]\n", buffer, infile.Name());
      return 1;
    }
    // Read number and freq
    if (infile.IO->Gets(buffer, BUFSIZE_)!=0) {
      mprinterr("Error: ReadEvecFile(): error while reading number and freq (%s)\n",infile.Name());
      return 1;
    }
    if (sscanf(buffer, "%i%lf", &nno, tmpval) != 2) {
      mprinterr("Error: ReadEvecFile(): error while scanning number and freq (%s)\n",infile.Name());
      return 1;
    }
    if (nno >= ibeg && nno <= iend)
      freq_[nvect_] = tmpval[0];
    else if (nno > iend)
      break;
    // Read coords
    nent = 0;
    for (int i = 0; i < nlines; ++i) {
      if (infile.IO->Gets(buffer, BUFSIZE_)!=0) {
        mprinterr("Error: ReadEvecFile(): error while reading evec coords (%s)\n",infile.Name());
        return 1;
      }
      nvals = sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf", tmpval, tmpval+1, tmpval+2,
                     tmpval+3, tmpval+4, tmpval+5, tmpval+6);
      for (int j = 0; j < nvals; j++) {
        evec_[ncoords * nvect_ + nent] = tmpval[j];
        ++nent;
      }
    }
    ++nvect_;
  }

  if (nvect_ != (iend - ibeg + 1)) {
    mprintf("Warning: Number of read evecs is %i, number of requested evecs is %i\n",
            nvect_, iend - ibeg + 1);
  }
 
  infile.CloseFile(); 
  return 0;
}

