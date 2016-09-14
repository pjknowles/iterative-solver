#ifndef PARAMETERVECTOR_H
#define PARAMETERVECTOR_H

#include <cstddef>
#include <vector>


namespace IterativeSolver {
  typedef double Scalar;

  class ParameterVector
  {
  public:
    ParameterVector(int variance=0);
    ParameterVector(Scalar* buffer, size_t length, int variance=0);
    virtual ~ParameterVector();
  private:
    int variance_;
    size_t length_;
    Scalar* buffer_;
  };

  class ParameterVectorSet : private std::vector<ParameterVector>
  {

  };

}

#endif // PARAMETERVECTOR_H
