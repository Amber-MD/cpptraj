// ProgressBar
#include "ProgressBar.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
// Takes value needed to be considered complete as input.
ProgressBar::ProgressBar(int maxIn) {
  max = maxIn;
  targetPercent=0;
}

// DESTRUCTOR
ProgressBar::~ProgressBar() { }

/* ProgressBar::Update()
 * If current percent is greater than target percent, print current
 * percent complete.
 */
void ProgressBar::Update(int current) {
  int currentPercent;

  currentPercent = (current * 100) / max;
  if (currentPercent >= targetPercent) {
    targetPercent+=10;
    mprintf("%2i%% ",currentPercent);
    mflush();
  }
}

/* ProgressBar::Complete()
 * Print 100%
 */
void ProgressBar::Complete() {
  mprintf("100%\n");
}

