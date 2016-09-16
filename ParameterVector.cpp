#include "ParameterVector.h"
#include <iostream>

using namespace IterativeSolver;

ParameterVector::ParameterVector(size_t length) :
    m_buffer(length), m_variance(0)
{
std::cout <<"ParameterVector bare constructor"<<std::endl;
}

ParameterVector::~ParameterVector() {}


ParameterVector& ParameterVector::operator=(const ParameterVector& other)
{
    m_buffer.resize(other.size());
    for (size_t k=0; k<size(); k++) m_buffer[k] = other[k];
    return *this;
}

void ParameterVector::axpy(ParameterScalar a, const ParameterVector& other)
{
    if (this->m_variance != other.m_variance) throw std::logic_error("mismatching co/contravariance");
    for (size_t k=0; k<size(); k++) m_buffer[k] += a*other[k];
}

void ParameterVector::zero()
{
    for (size_t k=0; k<size(); k++) m_buffer[k] = (ParameterScalar)0;
}

ParameterScalar ParameterVector::operator*(const ParameterVector& other) const
{
    if (this->m_variance * other.m_variance > 0) throw std::logic_error("mismatching co/contravariance");
    ParameterScalar result=0;
    for (size_t k=0; k<size(); k++) result += m_buffer[k] * other[k];
    return result;
}

