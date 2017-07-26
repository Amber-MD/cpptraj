#ifndef _xdr_seek_h
#define _xdr_seek_h

#ifdef __cplusplus
extern "C"
{
#endif

// for int64_t on older M$ Visual Studio
#if _MSC_VER && _MSC_VER < 1600 && !__INTEL_COMPILER
	#include "ms_stdint.h"
#else
	#include <stdint.h>
#endif

#include "xdrfile.h"

int64_t xdr_tell(XDRFILE *xd);
int xdr_seek(XDRFILE *xd, int64_t pos, int whence);

#ifdef __cplusplus
}
#endif

#endif
