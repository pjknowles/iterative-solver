#include "ParameterVector.h"
#include <iostream>

using namespace IterativeSolver;

ParameterVector::ParameterVector(int variance) :
    buffer_(nullptr), length_(0), variance_(variance)
{
std::cout <<"ParameterVector bare constructor"<<std::endl;
}
ParameterVector::ParameterVector(ParameterScalar* buffer, size_t length, int variance) :
    buffer_(buffer), length_(length), variance_(variance)
{
// std::cout << "ParameterVector constructor taking "<<*buffer<<std::endl;
// std::cout << "ParameterVector constructor set "<<*buffer_<<std::endl;
// std::cout << "ParameterVector constructor length_ "<<length_<<std::endl;
// std::cout << "ParameterVector constructor size() "<<size()<<std::endl;
}

ParameterVector::~ParameterVector() {}


ParameterVector& ParameterVector::operator=(const ParameterVector& other)
{
    length_=other.length_;
    buff.resize(length_);
    buffer_ = &buff[0];
    if (this->buffer_!=nullptr && other.buffer_!=nullptr) { // both in memory
        for (size_t k=0; k<length_; k++) buffer_[k] = other.buffer_[k];
//        std::cout << "ParameterVector operator= set "<<*buffer_<<std::endl;
    }
    else throw std::logic_error("implementation incomplete");
    std::cout << "ParameterVector::operator= has copied from "<<other.buffer_<<" to "<<buffer_<<std::endl;
    return *this;
}

void ParameterVector::axpy(ParameterScalar a, const ParameterVector& other)
{
    if (this->variance_ != other.variance_) throw std::logic_error("mismatching co/contravariance");
    if (this->buffer_!=nullptr && other.buffer_!=nullptr) { // both in memory
        for (size_t k=0; k<length_; k++) buffer_[k] += a*other.buffer_[k];
    }
    else throw std::logic_error("implementation incomplete");
}

void ParameterVector::zero()
{
  if (buffer_!=nullptr)
	for (size_t k=0; k<length_; k++) buffer_[k] = (ParameterScalar)0;
  else
        for (size_t k=0; k<length_; k++) buff[k] = (ParameterScalar)0;
}

ParameterScalar ParameterVector::operator*(const ParameterVector& other) const
{
    if (this->variance_ * other.variance_ > 0) throw std::logic_error("mismatching co/contravariance");
    ParameterScalar result=0;
    if (this->buffer_!=nullptr && other.buffer_!=nullptr) { // both in memory
        for (size_t k=0; k<length_; k++) std::cout << "in *, buffer elements:  "<< buffer_[k]<<" "<<other.buffer_[k]<<std::endl;
        for (size_t k=0; k<length_; k++) result += buffer_[k]*other.buffer_[k];
    }
    else throw std::logic_error("implementation incomplete");
    return result;
}

