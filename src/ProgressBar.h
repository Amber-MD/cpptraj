#ifndef INC_PROGRESSBAR_H
#define INC_PROGRESSBAR_H
/// Used to print progress to screen
class ProgressBar {
  public:
    ProgressBar();
    ProgressBar(int);

    void Update(int);
  private:
    static const int UNKNOWN_FRAMESIZE;
    int max_;
    float C_over_max_;
    float targetPercent_;
    bool unknownframes_;
};
#endif
