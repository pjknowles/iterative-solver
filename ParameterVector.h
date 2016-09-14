#ifndef PARAMETERVECTOR_H
#define PARAMETERVECTOR_H

#include <cstddef>
#include <vector>


namespace IterativeSolver {
  typedef double ParameterScalar;

  class ParameterVector
  {
  public:
    ParameterVector(int variance=0);
    ParameterVector(ParameterScalar* buffer, size_t length, int variance=0);
    virtual ~ParameterVector();
  private:
    int variance_;
    size_t length_;
    ParameterScalar* buffer_;
  };

  class ParameterVectorSet : private std::vector<ParameterVector>
  {

  };

}

#endif // PARAMETERVECTOR_H
