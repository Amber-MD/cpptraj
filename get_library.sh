#!/bin/bash
# Companion script to CPPTRAJ 'configure' script.
# Responsible for downloading, unpacking, configuring, and compiling external libraries.
# Daniel R. Roe
# 2021-03-03

WORKDIR=`pwd`

# Attempt to download and install a copy of library 
if [ -z "$SRCTAR" -o -z "$URL" -o -z "$LIBNAME" ] ; then
  echo "Error: Script variables are empty."
  exit 1
fi

# Check if --rebuild specified
REBUILD=0
if [ "$1" == '--rebuild' ] ; then
  REBUILD=1
  shift
fi

# First ask if we want to get the library
if [ $REBUILD -eq 0 ] ; then
  echo -n "Should CPPTRAJ attempt to build its own $LIBNAME? {y|n}: "
  read yesno
  while [ "$yesno" != 'y' -a "$yesno" != 'n' ] ; do
    echo -n "    Please enter 'y' or 'n': "
    read yesno
  done
  if [ "$yesno" = 'n' ] ; then
    exit 1
  fi
fi

# Get configure options
CONFIGOPTS=''
while [ ! -z "$1" ] ; do
  CONFIGOPTS="$CONFIGOPTS $1"
  shift
done

#echo "    CONFIGOPTS: $CONFIGOPTS"

CONFIGURE_LOG=$WORKDIR/$LIBNAME"_config.log"
COMPILE_LOG=$WORKDIR/$LIBNAME"_compile.log"

# Get library if necessary
if [ ! -f "$SRCTAR" ] ; then
  echo "    Downloading $LIBNAME..."
  WGET=`which wget`

  if [ -z "$WGET" ] ; then
    echo "Error: 'wget' not found. Cannot download $LIBNAME"
    exit 1
  fi

  $WGET $URL
  if [ $? -ne 0 -o ! -f "$SRCTAR" ] ; then
    echo "Error: Could not download $URL"
    exit 1
  fi
fi

# Get SRCDIR if necessary
if [ -z "$SRCDIR" ] ; then
  FIRSTFILE=`tar -tzf $SRCTAR | head -1`
  if [ ! -z "$FIRSTFILE" ] ; then
    SRCDIR=`dirname $FIRSTFILE`
    if [ "$SRCDIR" = '.' ] ; then
      SRCDIR=$FIRSTFILE
    fi
  fi
fi

# Unpack
if [ -z "$SRCDIR" ] ; then
  echo "Error: SRCDIR is empty."
  exit 1
fi
if [ ! -d "$SRCDIR" ] ; then
  echo "    Unpacking $LIBNAME..."
  tar -zxf $SRCTAR
  if [ $? -ne 0 -o ! -d "$SRCDIR" ] ; then
    echo "Error: Could not unpack $SRCTAR"
    exit 1
  fi
fi
if [ "$SRCDIR" = '.' -o ! -d "$SRCDIR" ] ; then
  echo "Error: $SRCDIR is not a directory."
  exit 1
fi
cd $SRCDIR

# Compiler/library-specific modifications
MAKE_TARGET=''
if [ -f 'make.inc' ] ; then
  rm 'make.inc'
fi
if [ "$LIBNAME" = 'lapack' ] ; then
  MAKE_TARGET='blaslib lapacklib'
  lapackflags=''
  if [ "$FC" = 'gfortran' ] ; then
    lapackflags='-frecursive'
  fi
  # LAPACK has no configure; requires make.inc
  cat > make.inc <<EOF
SHELL = /bin/sh
FORTRAN  = $FC
OPTS     = -O2 $lapackflags
DRVOPTS  = \$(OPTS)
NOOPT    = -O0 $lapackflags
LOADER   = $FC
LOADOPTS =
#TIMER    = INT_ETIME
TIMER    = INT_CPU_TIME
CC           = $CC
CFLAGS       = $CFLAGS 
ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib
XBLASLIB     =
BLASLIB      = libblas.a
CBLASLIB     = 
LAPACKLIB    = liblapack.a
TMGLIB       = libtmglib.a
LAPACKELIB   = liblapacke.a
EOF
fi

# Configure
#echo ""
#echo "    CC=$CC"
#echo "    CFLAGS=$CFLAGS"
#echo "    PREFIX=$PREFIX"
#echo "    FC=$FC"
#echo "    FFLAGS=$FFLAGS"
echo -n "    Configuring $LIBNAME... "
if [ "$LIBNAME" = 'lapack' ] ; then
  echo -n "(using generated make.inc) "
elif [ -f 'configure' ] ; then
  # Run configure
  CC="$CC" CFLAGS="$CFLAGS" FC="$FC" FFLAGS="$FFLAGS" ./configure --prefix=$PREFIX $CONFIGOPTS > $CONFIGURE_LOG 2>&1
  if [ $? -ne 0 ] ; then
    echo "Failed."
    echo "Check $CONFIGURE_LOG for errors."
    exit 1
  fi
elif [ -f 'Makefile' ] ; then
  # No configure - try to modify Makefile on the fly.
  if [ ! -f 'Makefile.original' ] ; then
    cp Makefile Makefile.original
  fi
  awk -v cc="$CC" -v cflags="$CFLAGS" -v fc="$FC" -v fflags="$FFLAGS" -v prefix="$PREFIX" '{
    if (index($1,"CC=")!=0)
      printf("CC=%s\n", cc);
    else if (index($1,"FC=")!=0)
      printf("FC=%s\n", fc);
    else if (index($1,"CFLAGS=")!=0)
      printf("CFLAGS=%s\n", cflags);
    else if (index($1,"FFLAGS=")!=0)
      printf("FFLAGS=%s\n", fflags);
    else if (index($1,"PREFIX=")!=0)
      printf("PREFIX=%s\n", prefix);
    else
      print $0;
  }' Makefile.original > Makefile
  if [ $? -ne 0 ] ; then
    exit 1
  fi
else
  echo "Error: $LIBNAME has no 'configure' or 'Makefile'."
  exit 1
fi
echo "Success."

# Determine make command
if [ -z "$MAKE_COMMAND" ] ; then
  NPROC=`nproc`
  if [ -z "$NPROC" ] ; then
    MAKE_COMMAND='make'
  elif [ $NPROC -le 2 ] ; then
    MAKE_COMMAND='make'
  else
    HALF=`echo "$NPROC / 2" | bc`
    MAKE_COMMAND="make -j$HALF"
  fi
  echo "    MAKE_COMMAND is not set; set to '$MAKE_COMMAND'"
fi

# Build
echo -n "    Compiling $LIBNAME (may be time-consuming)... "
make clean > $COMPILE_LOG 2>&1
$MAKE_COMMAND $MAKE_TARGET > $COMPILE_LOG 2>&1
if [ $? -ne 0 ] ; then
  echo "Build failed."
  echo "Check $COMPILE_LOG for errors."
  exit 1
fi

# Install
if [ "$LIBNAME" = 'lapack' ] ; then
  # Only made blas and lapack libraries; need to move them manually
  if [ -f 'BLAS/SRC/libblas.a' ] ; then
    blaslib=BLAS/SRC/libblas.a
  else
    echo "Error: BLAS not made."
    exit 1
  fi
  if [ -f 'liblapack.a' ] ; then
    lapacklib=liblapack.a
  elif [ -f 'SRC/liblapack.a' ] ; then
    lapacklib=SRC/liblapack.a
  else
    echo "Error: LAPACK not made."
    exit 1
  fi
  mv $blaslib $lapacklib $PREFIX/lib/
else
  make install >> $COMPILE_LOG 2>&1
fi
if [ $? -ne 0 ] ; then
  echo "Install failed."
  echo "Check $COMPILE_LOG for errors."
  exit 1
fi

echo "Success."

exit 0
