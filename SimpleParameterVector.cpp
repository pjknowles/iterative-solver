#include "SimpleParameterVector.h"

#include <stdexcept>
#include <iostream>

using namespace IterativeSolver;

SimpleParameterVector::SimpleParameterVector(size_t length) :
    m_buffer(length), m_variance(0)
{
    std::cout << "SimpleParameterVector constructed, size()="<<size()<<std::endl;
}

SimpleParameterVector::~SimpleParameterVector() {}


SimpleParameterVector& SimpleParameterVector::operator=(const SimpleParameterVector& other)
{
    m_buffer.resize(other.m_buffer.size());
    for (size_t k=0; k<m_buffer.size(); k++) m_buffer[k] = other.m_buffer[k];
    return *this;
}

void SimpleParameterVector::axpy(ParameterScalar a, const SimpleParameterVector& other)
{
    if (this->m_variance != other.m_variance) throw std::logic_error("mismatching co/contravariance");
    for (size_t k=0; k<m_buffer.size(); k++) m_buffer[k] += a*other.m_buffer[k];
}

void SimpleParameterVector::zero()
{
    for (size_t k=0; k<m_buffer.size(); k++) m_buffer[k] = (ParameterScalar)0;
}

ParameterScalar SimpleParameterVector::operator*(const SimpleParameterVector& other) const
{
    if (this->m_variance * other.m_variance > 0) throw std::logic_error("mismatching co/contravariance");
    ParameterScalar result=0;
    for (size_t k=0; k<m_buffer.size(); k++) result += m_buffer[k] * other.m_buffer[k];
    return result;
}


  ParameterScalar& SimpleParameterVector::operator[](size_t pos) { return this->m_buffer[pos];}
  const ParameterScalar& SimpleParameterVector::operator[](size_t pos) const { return this->m_buffer[pos];}
  size_t SimpleParameterVector::size() const { return this->m_buffer.size();}
