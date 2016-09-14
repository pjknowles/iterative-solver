#include "ParameterVector.h"

using namespace IterativeSolver;

ParameterVector::ParameterVector(int variance) :
  buffer_(nullptr), length_(0), variance_(variance)
{

}
ParameterVector::ParameterVector(ParameterScalar* buffer, size_t length, int variance) :
  buffer_(buffer), length_(length), variance_(variance)
{

}

ParameterVector::~ParameterVector() {}


ParameterVector& ParameterVector::operator=(const ParameterVector& other)
{
 return *this;
}

ParameterVector& ParameterVector::operator *=(ParameterScalar a)
{
    if (this->buffer_!=nullptr) { //both in memory
        for (size_t k=0; k<length_; k++) buffer_[k] *= a;
    }
    else throw std::logic_error("implementation incomplete");
 return *this;
}

ParameterScalar ParameterVector::operator*(const ParameterVector& other) const
{
    ParameterScalar result=0;
    if (this->buffer_!=nullptr && other.buffer_!=nullptr) { // both in memory
        for (size_t k=0; k<length_; k++) result += buffer_[k]*other.buffer_[k];
    }
    else throw std::logic_error("implementation incomplete");
    return result;
}

