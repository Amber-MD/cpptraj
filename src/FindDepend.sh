#!/bin/bash

# For each cpp file list the include files it depends on
# Create an entry like: <obj file> : <cpp file> <include file1> ...
for FILE in `ls *.cpp *.c *.h` ; do
  CPPFILE=$FILE
  FILE=${FILE/%cpp/o}
  FILE=${FILE/%c/o}
  #echo -n "$FILE : $CPPFILE"
  awk -v ofile=$FILE -v cppfile=$CPPFILE 'BEGIN{
    Ndepend=0;
    Dependline="";
    if (ofile==cppfile)
      Dependline=ofile" :";
    else
      Dependline=ofile" : "cppfile;
  }{
    if ($1=="#include") 
      if ( index($2,"<")==0 && $2!="\"mpi.h\"" && $2!="\"zlib.h\"" && $2!="\"bzlib.h\"" && $2!="\"netcdf.h\"" ) {
        gsub(/"/,"",$2);
        #printf(" %s",$2);
        Dependline=Dependline" "$2;
        Ndepend++;
      }
}END{
  if (Ndepend>0) print Dependline;
  #printf("\n");
}' $CPPFILE
done

