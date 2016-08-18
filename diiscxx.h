#ifndef DIISCXX_H
#define DIISCXX_H
#include <string>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>


class diiscxx
{
public:
  enum DiisMode_type {disabled, DIIS, KAIN};
  diiscxx(size_t AmpLength, size_t ResLength=0, size_t OtherLength=0, double threshold=1e-3, DiisMode_type DiisMode=DIIS);
  virtual void dumpInit();
  virtual void dump(const double* buffer, size_t length, size_t address);
  virtual void load(double* buffer, int index, size_t length, size_t address);
  virtual void dumpEnd();
  void setDiisMode(enum DiisMode_type DiisMode =DIIS);
private:
  diiscxx();
  size_t AmpLength_;
  size_t ResLength_;
  size_t OtherLength_;
  size_t TotalLength_;
  enum DiisMode_type DiisMode_;
  double threshold_;

  // part of I/O implementation
  std::fstream dumpFile;
};

#endif // DIISCXX_H
