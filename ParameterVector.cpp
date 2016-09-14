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
    length_=other.length_;
    buff.resize(length_);
    buffer_ = &buff[0];
    if (this->buffer_!=nullptr && other.buffer_!=nullptr) { // both in memory
        for (size_t k=0; k<length_; k++) buffer_[k] = other.buffer_[k];
    }
    else throw std::logic_error("implementation incomplete");
    return *this;
}

void ParameterVector::axpy(ParameterScalar a, ParameterVector& other)
{
    if (this->variance_ != other.variance_) throw std::logic_error("mismatching co/contravariance");
    if (this->buffer_!=nullptr && other.buffer_!=nullptr) { // both in memory
        for (size_t k=0; k<length_; k++) buffer_[k] += a*other.buffer_[k];
    }
    else throw std::logic_error("implementation incomplete");
}

ParameterScalar ParameterVector::operator*(const ParameterVector& other) const
{
    if (this->variance_ * other.variance_ > 0) throw std::logic_error("mismatching co/contravariance");
    ParameterScalar result=0;
    if (this->buffer_!=nullptr && other.buffer_!=nullptr) { // both in memory
        for (size_t k=0; k<length_; k++) result += buffer_[k]*other.buffer_[k];
    }
    else throw std::logic_error("implementation incomplete");
    return result;
}

