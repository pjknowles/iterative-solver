#include "ParameterVector.h"
#include <stdexcept>
#include <iostream>

using namespace IterativeSolver;

ParameterVector::ParameterVector(size_t length) :
    m_buffer(length), m_variance(0)
{
}

ParameterVector::ParameterVector(const ParameterVector &source)
{
//    std::cout << "ParameterVector:: copy constructor"<<std::endl;
    *this = source;
//    this->operator=(source);
//    std::cout << "ParameterVector:: copy constructor, size()="<<size()<<std::endl;
}

ParameterVector::~ParameterVector() {}


ParameterVector& ParameterVector::operator=(const ParameterVector& other)
{
    m_buffer.resize(other.m_buffer.size());
    for (size_t k=0; k<m_buffer.size(); k++) m_buffer[k] = other.m_buffer[k];
    m_variance = other.m_variance;
//    std::cout << "ParameterVector::operator=, size()="<<size()<<std::endl;
    return *this;
}

ParameterVector* ParameterVector::clone() const
{
//    std::cout <<"ParameterVector::clone *this="<<*this<<std::endl;
  ParameterVector* result = new ParameterVector(*this);
//    std::cout <<"ParameterVector::clone *result="<<*result<<std::endl;
  return new ParameterVector(*this);
}

void ParameterVector::axpy(ParameterScalar a, const ParameterVector *other)
{
    if (this->m_variance != other->m_variance) throw std::logic_error("mismatching co/contravariance");
    for (size_t k=0; k<m_buffer.size(); k++) m_buffer[k] += a*other->m_buffer[k];
}

void ParameterVector::zero()
{
    for (size_t k=0; k<m_buffer.size(); k++) m_buffer[k] = (ParameterScalar)0;
}

ParameterScalar ParameterVector::dot(const ParameterVector* other) const
{
    if (this->m_variance * other->m_variance > 0) throw std::logic_error("mismatching co/contravariance");
    ParameterScalar result=0;
//    std::cout << "ParameterVector::dot"<<std::endl;
    for (size_t k=0; k<m_buffer.size(); k++) result += m_buffer[k] * other->m_buffer[k];
    return result;
}



  ParameterScalar& ParameterVector::operator[](size_t pos) { return this->m_buffer[pos];}
  const ParameterScalar& ParameterVector::operator[](size_t pos) const { return this->m_buffer[pos];}
  size_t ParameterVector::size() const { return this->m_buffer.size();}



