#ifndef INC_PROGRESSBAR_H
#define INC_PROGRESSBAR_H

class ProgressBar {
    int max;
    int targetPercent;
  public:
    ProgressBar(int);
    ~ProgressBar();

    void Update(int);
    void Complete();
};

#endif
