// ProgressBar
#include "ProgressBar.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
/// Takes value needed to be considered complete as input.
ProgressBar::ProgressBar(int maxIn) {
  max = maxIn - 1;
  oneframe=false;
  unknownframes=false;
  C_over_max = 1;
  targetPercent=0;
  if (max == 0)
    oneframe=true;
  else if (max < 0) {
    unknownframes=true;
    mprintf("Progress: '+' = 200 frames processed.\n");
    targetPercent=199;
  } else {
    C_over_max = (float) max;
    C_over_max = 100 / C_over_max;
  }
  first=true;
}

// ProgressBar::Update()
/** If current percent is greater than target percent, print current
  * percent complete.
  */
void ProgressBar::Update(int current) {
  float currentPercent;
  int target;

  if (unknownframes) {
    if (first) {
      mprintf("%10i ",current);
      first=false;
      mflush();
    }
    currentPercent = (float) current;
    if (currentPercent > targetPercent) {
      mprintf("+");
      target = (int) targetPercent;
      target++;
      if ((target % 5000) == 0)
        mprintf("\n%10i ",current);
      targetPercent+=200;
      mflush();
    } 
  } else {
    currentPercent = current * C_over_max;
    if (currentPercent >= targetPercent) {
      if (current==max) 
        mprintf("100%% Complete.\n");
      else {
        targetPercent+=10;
        mprintf("%2.0f%% ",currentPercent);
      }
      //if (first) {
      //  mflush();
      //  first=false;
      //}
      mflush();
    }
  }
}

// ProgressBar::PrintBar()
/** Print an actual progress bar based on current frame. */
void ProgressBar::PrintBar(int current) {
  char buffer[128];
  int i,j,percent,numChars;

  if (unknownframes) return; //NOTE: Eventually print dots like ptraj?

  // Fraction complete, max 100. In case of 1 frame (stop=0) set 100
  if (oneframe)
    percent = 100;
  else
    percent = (current * 100) / max;
    //percent = current * C_over_max;

  buffer[0]='[';
  i=1;
  // Set number of characters
  numChars = percent / 2;
  //mprintf("[%i] %i : numChars=%i\n",worldrank,set,numChars);
  for (j=0; j<numChars; j++) buffer[i++]='|';
  // Fill the rest
  for (; j<50; j++)
    buffer[i++]='-';
  // Add last character and reset/newline character
  if (current==max) buffer[i-1]='|';
  buffer[i++]=']';
  buffer[i++]='\r';
  if (current==max) {
    buffer[i-1]=' ';
    buffer[i++]='C';
    buffer[i++]='o';
    buffer[i++]='m';
    buffer[i++]='p';
    buffer[i++]='l';
    buffer[i++]='e';
    buffer[i++]='t';
    buffer[i++]='e';
    buffer[i++]='.';
    buffer[i++]='\n';
  }
  // Finish off and print
  buffer[i]='\0';
  mprintf("  %s",buffer);
  // On first frame flush so that progress bar appears
  if (first) { mflush(); first=false;}

/*
  // DEBUG - print without the restore char
  fprintf(stderr,"ProgressBar: stop=%i currentFrame=%i percent=%i\n",
          stop,currentFrame,percent);
  for (j=0; j < i; j++)
    if (buffer[j]=='\r') buffer[j]='\n';
  fprintf(stderr,"           : %s\n",buffer);
*/
}

