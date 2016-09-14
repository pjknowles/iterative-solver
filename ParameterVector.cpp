#include "ParameterVector.h"

using namespace IterativeSolver;

ParameterVector::ParameterVector(int variance) :
  buffer_(nullptr), length_(0), variance_(variance)
{

}
ParameterVector::ParameterVector(Scalar* buffer, size_t length, int variance) :
  buffer_(buffer), length_(length), variance_(variance)
{

}

ParameterVector::~ParameterVector() {}

