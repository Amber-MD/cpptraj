#include "NDiff.h"
#include "ArgList.h"
#include "BufferedLine.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"
#include <cmath>
#include <algorithm>

/// Return the maximum relative error of two values
double maxrelerr(double x, double y)
{
    if (x == y)
	return (0);
    else if ((x != 0) && (y != 0))
	return (fabs(x-y) / std::min(fabs(x), fabs(y)));
    else if ((x == 0) && (y != 0))
	return (1);
    else if ((y == 0) && (x != 0))
	return (1);
    else
	return (0);
}

void report_difference(CpptrajFile& outfile, int& ndiff, int NRLINE,
                       const char* f1line, const char* f2line)
{
  outfile.Printf("%dc%d\n", NRLINE, NRLINE);
  outfile.Printf("< %s\n", f1line);
  outfile.Printf("> %s\n", f2line);
  ndiff++;
}

static double ndiff_max_abserr = 0;
static int    ndiff_max_abserr_line = 0;
static int    ndiff_max_abserr_field = 0;
static double ndiff_max_relerr = 0;
static int    ndiff_max_relerr_line = 0;
static int    ndiff_max_relerr_field = 0;

int diff_field(std::string const& f1, std::string const& f2, bool is_absolute, double tolIn,
               int NRLINE, int nfield)
{
   // If both fields are identical as strings, return 0.
   if (f1 == f2) return 0;
   // If both fields are numeric, return 0 if they are close enough 
   if (validDouble(f1) && validDouble(f2)) {
     double d1 = convertToDouble(f1);
     double d2 = convertToDouble(f2);
     double This_Abserr = fabs(d1 - d2);
     double This_Relerr = maxrelerr(d1, d2);

     bool is_diff = false;
     if (is_absolute) {
       if ( This_Abserr > tolIn )
         is_diff = true;
     } else {
       if ( This_Relerr > tolIn )
         is_diff = true;
     }
     if (is_diff) {
       if (This_Abserr > ndiff_max_abserr) {
         ndiff_max_abserr = This_Abserr;
         ndiff_max_abserr_line = NRLINE;
         ndiff_max_abserr_field = nfield;
       }
       if (This_Relerr > ndiff_max_relerr) {
         ndiff_max_relerr = This_Relerr;
         ndiff_max_relerr_line = NRLINE;
         ndiff_max_relerr_field = nfield;
       }
       return 1;
     } else
       return 0;
  } else
    return 1;
}

int compare_all( CpptrajFile& outfile, int& ndiff, int NRLINE,
                 ArgList const& f1line, ArgList const& f2line,
                 bool is_absolute, double tolIn )
{
  for (int iarg = 0; iarg < f1line.Nargs(); iarg++) {
    if (diff_field(f1line[iarg], f2line[iarg], is_absolute, tolIn, NRLINE, iarg)) {
      report_difference(outfile, ndiff, NRLINE, f1line.ArgLine(), f2line.ArgLine());
      return 1;
    }
  }
  return 0;
}

int ndiff_compare_files(BufferedLine& file1, BufferedLine& file2, bool is_absolute,
                        double tolIn)
{
  int NRLINE = 0;
  int ndiff = 0;

  CpptrajFile outfile;
  outfile.OpenWrite("");

  const char* ptr1 = file1.Line();
  const char* ptr2 = file2.Line();

  while (ptr1 != 0 && ptr2 != 0) {
    NRLINE++;
    ArgList f1line( ptr1 );
    ArgList f2line( ptr2 );
    if (f1line.Nargs() == f2line.Nargs()) {
      compare_all( outfile, ndiff, NRLINE, f1line, f2line, is_absolute, tolIn );
    } else {
      report_difference(outfile, ndiff, NRLINE, ptr1, ptr2 );
    }
    ptr1 = file1.Line();
    ptr2 = file2.Line();
  }

  //if (QUIET == 0)
  //{
    if (ndiff_max_abserr > 0)
      outfile.Printf("### Maximum absolute error in matching lines = %.2e at line %d field %d\n",
                     ndiff_max_abserr, ndiff_max_abserr_line, ndiff_max_abserr_field+1);
    if (ndiff_max_relerr > 0)
      outfile.Printf("### Maximum relative error in matching lines = %.2e at line %d field %d\n",
                     ndiff_max_relerr, ndiff_max_relerr_line, ndiff_max_relerr_field+1);
  //}
  if (ptr2 != 0) {
    mprintf("Warning: file %s is short.\n", file1.Filename().full());
    ndiff++;
  }
  if (ptr1 != 0) {
    mprintf("Warning: file %s is short.\n", file2.Filename().full());
    ndiff++;
  }

  return ndiff;
}

/** Intended to be a faster drop-in replacement for Nelson H. F. Beebe's 
  * ndiff.awk script.
  */
int NDiff(std::string const& tolArgStr, std::string const& fname1, std::string const& fname2)
{
  // Parse the tolerance arg
  ArgList tolarg(tolArgStr, "=");
  if (tolarg.Nargs() != 2) {
    mprinterr("Error: ndiff: malformed tolerance arg: %s\n", tolArgStr.c_str());
    return 1;
  }
  
  bool is_absolute = true;
  if (tolarg[0] == "RELERR")
    is_absolute = false;
  else if (tolarg[0] == "ABSERR")
    is_absolute = true;
  else {
    mprinterr("Error: ndiff: Expected 'RELERR' or 'ABSERR', got '%s'\n", tolarg[0].c_str());
    return -1;
  }

  if (!validDouble(tolarg[1])) {
    mprinterr("Error: ndiff: '%s' is not a valid tolerance.\n", tolarg[1].c_str());
    return -1;
  }
  double tolIn = convertToDouble(tolarg[1]);

  BufferedLine file1, file2;
  if (file1.OpenFileRead( fname1 )) {
    mprinterr("Error: ndiff: Could not open '%s'\n", fname1.c_str());
    return -1;
  }
  if (file2.OpenFileRead( fname2 )) {
    mprinterr("Error: ndiff: Could not open '%s'\n", fname2.c_str());
    return -1;
  }

  return ndiff_compare_files(file1, file2, is_absolute, tolIn);
}
