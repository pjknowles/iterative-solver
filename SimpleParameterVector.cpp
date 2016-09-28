#include "SimpleParameterVector.h"

#include <stdexcept>
#include <iostream>

using namespace IterativeSolver;

SimpleParameterVector::SimpleParameterVector(size_t length) :
  m_buffer(length), m_variance(0)
{
}

SimpleParameterVector::SimpleParameterVector(const SimpleParameterVector &source)
{
  *this = source;
}

SimpleParameterVector::~SimpleParameterVector() {}


SimpleParameterVector& SimpleParameterVector::operator=(const SimpleParameterVector& other)
{
  m_buffer.resize(other.m_buffer.size());
  for (size_t k=0; k<m_buffer.size(); k++) m_buffer[k] = other.m_buffer[k];
  m_variance = other.m_variance;
  return *this;
}

SimpleParameterVector* SimpleParameterVector::clone() const
{
  SimpleParameterVector* result = new SimpleParameterVector(*this);
  return new SimpleParameterVector(*this);
}

void SimpleParameterVector::axpy(ParameterScalar a, const ParameterVector* other)
{
  const SimpleParameterVector* othe=dynamic_cast <const SimpleParameterVector*> (other);
  if (this->m_variance != othe->m_variance) throw std::logic_error("mismatching co/contravariance");
  for (size_t k=0; k<m_buffer.size(); k++) m_buffer[k] += a*othe->m_buffer[k];
}

void SimpleParameterVector::zero()
{
  for (size_t k=0; k<m_buffer.size(); k++) m_buffer[k] = (ParameterScalar)0;
}

ParameterScalar SimpleParameterVector::dot(const ParameterVector *other) const
{
  const SimpleParameterVector* othe=dynamic_cast <const SimpleParameterVector*> (other);
  if (this->m_variance * othe->m_variance > 0) throw std::logic_error("mismatching co/contravariance");
  ParameterScalar result=0;
  for (size_t k=0; k<m_buffer.size(); k++) result += m_buffer[k] * othe->m_buffer[k];
  return result;
}


ParameterScalar& SimpleParameterVector::operator[](size_t pos) { return this->m_buffer[pos];}
const ParameterScalar& SimpleParameterVector::operator[](size_t pos) const { return this->m_buffer[pos];}
size_t SimpleParameterVector::size() const { return this->m_buffer.size();}
