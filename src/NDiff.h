#ifndef INC_NDIFF_H
#define INC_NDIFF_H
#include <string>
/// Tolerance arg (ABSERR=<tol> or RELERR=<tol>), file1, file2
int NDiff(std::string const&, std::string const&, std::string const&);

#endif
