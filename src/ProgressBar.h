#ifndef INC_PROGRESSBAR_H
#define INC_PROGRESSBAR_H
/// Used to print progress to screen
class ProgressBar {
  public:
    ProgressBar(int);

    void Update(int);
    void PrintBar(int);
  private:
    int max;
    float C_over_max;
    float targetPercent;
    bool first;
    bool oneframe;
    bool unknownframes;
};
#endif
