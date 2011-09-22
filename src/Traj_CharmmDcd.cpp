// CharmmDcd
#include "Traj_CharmmDcd.h"
#include "CpptrajStdio.h"
#include <cstddef>

// CONSTRUCTOR
CharmmDcd::CharmmDcd() {
  dcdatom=0;
}

// DESTRUCTOR
CharmmDcd::~CharmmDcd() {
}

/* CharmmDcd::openTraj()
 */
int CharmmDcd::openTraj() {
  int err;

  err=0;
  switch (tfile->access) {
    case READ : err = tfile->OpenFile(); break;
    case APPEND :
      mprintf("Error: Append not supported for dcd files.\n");
      err=1;
      break;
    case WRITE :
      // Set up title
      if (title==NULL)
        this->SetTitle((char*)"Cpptraj generated mol2 file.\0");
      err = tfile->OpenFile();
      break;
  }
      
  return err;
}

/* CharmmDcd::closeTraj() {
 */
void CharmmDcd::closeTraj() {
  tfile->CloseFile();
}

// Convert 4 bytes to integer
/*static int byte_to_int(unsigned char byte[4]) {
  int val,temp;
  val=0;
  temp = (int) byte[0];
  temp <<= 24;
  val = val | temp;
  temp = (int) byte[1];
  temp <<= 16;
  val = val | temp;
  temp = (int) byte[2];
  temp <<= 8;
  val = val | temp;
  temp = (int) byte[3];
  val = val | temp;
  return val;
}*/

int CharmmDcd::readBinaryInteger(int byteorder) {
  byte u;
  if (byteorder) {
    tfile->IO->Read(u.c+3,1,1);
    tfile->IO->Read(u.c+2,1,1);
    tfile->IO->Read(u.c+1,1,1);
    tfile->IO->Read(u.c  ,1,1);
  } else {
    tfile->IO->Read(u.c  ,1,1);
    tfile->IO->Read(u.c+1,1,1);
    tfile->IO->Read(u.c+2,1,1);
    tfile->IO->Read(u.c+3,1,1);
  }
  return(u.i);
}


/* CharmmDcd::binaryByteOrder()
 * ADAPTED FROM PTRAJ:
 *  Routine to read in 4 bytes from the FILE *fd and order in both
 *  big and little endian format...  If the value is "expected",
 *  the byteorder is set (1: little, 0: big).  If the value for both
 *  is not expected, the smaller of the two orders is returned (but a
 *  warning is given).
 */
int CharmmDcd::binaryByteOrder(int expected, int *byteorder) {
  int i_bigE,i_litE;
  //unsigned char byte_bigE[4];
  //unsigned char byte_litE[4];
  byte byte0;
  byte byte1;
 
  // Load as Big Endian
  //tfile->IO->Read(byte_bigE  ,1,1);
  //tfile->IO->Read(byte_bigE+1,1,1);
  //tfile->IO->Read(byte_bigE+2,1,1);
  //tfile->IO->Read(byte_bigE+3,1,1);
  //i_bigE = byte_to_int(byte_bigE);
  tfile->IO->Read(byte0.c  ,1,1);
  tfile->IO->Read(byte0.c+1,1,1);
  tfile->IO->Read(byte0.c+2,1,1);
  tfile->IO->Read(byte0.c+3,1,1);
  i_bigE = byte0.i;

  // Flip to Little Endian
  //byte_litE[3] = byte_bigE[0];
  //byte_litE[2] = byte_bigE[1];
  //byte_litE[1] = byte_bigE[2];
  //byte_litE[0] = byte_bigE[3];
  //i_litE = byte_to_int(byte_litE);
  byte1.c[3] = byte0.c[0];
  byte1.c[2] = byte0.c[1];
  byte1.c[1] = byte0.c[2];
  byte1.c[0] = byte0.c[3];
  i_litE = byte1.i;

  if (debug >= 0) {
    mprintf("SGI      (big endian) order value is %d\n", i_bigE);
    mprintf("Linux (little endian) order value is %d\n", i_litE);
  }

  if (i_bigE == expected) {
    *byteorder = 0;
    return i_bigE;
  } else if (i_litE == expected) {
    *byteorder = 1;
    return i_litE;
  }

  /*  If an unexpected value is found (i.e. not == expected), then don't die,
   *  just return the smaller of the two values and set the endian-ness as 
   *  appropriate...
   */
  mprintf("binaryByteOrder() [binary I/O]: integer value read was not the\n");
  mprintf("expected value of %i.  Big endian: %i  Little endian: %i\n",
          expected, i_bigE, i_litE);
  mprintf("Returning the smaller value, endian-ness\n");
  
  if (i_bigE < i_litE) {
    *byteorder = 0;
    return i_bigE;
  } 
  *byteorder = 1;
  return i_litE;
}

/* CharmmDcd::setupRead()
 */
int CharmmDcd::setupRead(AmberParm *trajParm) {
  int Frames=0;
  int byteorder=0;
  int icntrl[20];
  byte magic;

  if ( openTraj() ) return -1;

  binaryByteOrder(84,&byteorder);

  magic.i = readBinaryInteger(0);
  mprintf("CHARMMDCD: Magic %c %c %c %c\n",magic.c[0],magic.c[1],magic.c[2],magic.c[3]);

  // load ICNTRL variables and print for debugging
  for (int i=0; i < 20; i++) { 
    icntrl[i] = readBinaryInteger(byteorder);
    mprintf("\ticntrl[%i]= %i\n",i,icntrl[i]);
  }

  if (icntrl[11] == 1) {
    mprintf("Warning: CharmmDcd::setupRead(): DIM4 (fourth dimension) not supported yet!\n");
    return -1;
  }


  closeTraj();

  return Frames;
}

/* CharmmDcd::readFrame()
 */
int CharmmDcd::readFrame(int set,double *X, double *V,double *box, double *T) {

  return 0;
}

/* CharmmDcd::processWriteArgs()
 */
int CharmmDcd::processWriteArgs(ArgList *argIn) {
  return 0;
}

/* CharmmDcd::setupWrite()
 * Set parm information required for write, and check write mode against
 * number of frames to be written.
 */
int CharmmDcd::setupWrite(AmberParm *trajParm) {
  return 0;
}

/* CharmmDcd::writeFrame()
 */
int CharmmDcd::writeFrame(int set, double *X, double *V,double *box, double T) {

  return 0;
}
 
/* CharmmDcd::info()
 */
void CharmmDcd::info() {
  mprintf("is a CHARMM DCD file");
}
