# Source this from all unit test scripts
. ../../MasterTest.sh

# The variable UNITSOURCES will list all cpptraj source files needed by the unit test.

# Create Makefile
CreateMakefile() {
  cat > Makefile <<EOF
include $CPPTRAJ_CONFIGH

main.o : main.cpp
	\$(CXX) \$(CXXFLAGS) -I$CPPTRAJ_SRCDIR -c -o main.o main.cpp
EOF
  UNITOBJECTS='main.o'
  for CSRCFILE in $UNITSOURCES ; do
    # Create object name
    ONAME=${CSRCFILE/.cpp/.o}
    UNITOBJECTS="$UNITOBJECTS $ONAME"
    # Add the rule
  cat >> Makefile <<EOF

$ONAME: $CPPTRAJ_SRCDIR/$CSRCFILE
	\$(CXX) \$(CXXFLAGS) -I$CPPTRAJ_SRCDIR -c -o $ONAME $CPPTRAJ_SRCDIR/$CSRCFILE
EOF
  done
  cat >> Makefile <<EOF
all: a.out

a.out: $UNITOBJECTS
	\$(CXX) -o a.out $UNITOBJECTS

EOF
}

# Run make, check result
RunMake() {
  ((PROGCOUNT++))
  echo ""
  echo "  CPPTRAJ: $1"
  if [ -z "$CPPTRAJ_DACDIF" ] ; then
    OUT "  CPPTRAJ: $1"
  fi
  make all
  if [ $? -ne 0 -o ! -f 'a.out' ] ; then
    ProgramError "$1"
  else
    if [ -z "$CPPTRAJ_ERROR" ] ; then
      ./a.out
    else
      $VALGRIND ./a.out 2>> $CPPTRAJ_ERROR
    fi
    if [ $? -ne 0 ] ; then
      ProgramError "$1"
    fi
  fi
}

if [ -z "$CPPTRAJHOME" -o ! -d "$CPPTRAJHOME" ] ; then
  echo "Cannot find CPPTRAJ home directory, set CPPTRAJHOME."
  exit 1
fi

# Determine cpptraj source directory, CPPTRAJ_SRCDIR
if [ -z "$CPPTRAJ_SRCDIR" ] ; then

  CPPTRAJ_SRCDIR=$CPPTRAJHOME/src
  echo "  CPPTRAJ source directory: $CPPTRAJ_SRCDIR"
  if [ ! -d "$CPPTRAJ_SRCDIR" ] ; then
    echo "CPPTRAJ source directory not found."
    exit 1
  fi

  if [ ! -f "$CPPTRAJ_SRCDIR/main.cpp" ] ; then
    echo "CPPTRAJ source directory $CPPTRAJ_SRCDIR appears empty."
    exit 1
  fi
  export CPPTRAJ_SRCDIR
fi

# Determine the config.h file, CPPTRAJ_CONFIGH
if [ -z "$CPPTRAJ_CONFIGH" ] ; then
  CPPTRAJ_CONFIGH=$CPPTRAJHOME/config.h
  echo "  CPPTRAJ config.h file: $CPPTRAJ_CONFIGH"
  if [ ! -f "$CPPTRAJ_CONFIGH" ] ; then
    echo "CPPTRAJ config.h file not found."
    exit 1
  fi
  export CPPTRAJ_CONFIGH
fi

