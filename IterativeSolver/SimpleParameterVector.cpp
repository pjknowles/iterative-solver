#include "SimpleParameterVector.h"

#include <stdexcept>
#include <iostream>

using namespace IterativeSolver;

SimpleParameterVector::SimpleParameterVector(size_t length) :
  m_buffer(length)
{
}

SimpleParameterVector::~SimpleParameterVector() {}


SimpleParameterVector& SimpleParameterVector::operator=(const SimpleParameterVector& other)
{
  m_buffer.resize(other.m_buffer.size());
  for (size_t k=0; k<m_buffer.size(); k++) m_buffer[k] = other.m_buffer[k];
  setVariance(other.variance());
  return *this;
}

void SimpleParameterVector::axpy(ParameterScalar a, const ParameterVector* other)
{
  const SimpleParameterVector* othe=dynamic_cast <const SimpleParameterVector*> (other);
  if (this->variance() != othe->variance()) throw std::logic_error("mismatching co/contravariance");
  for (size_t k=0; k<m_buffer.size(); k++) m_buffer[k] += a*othe->m_buffer[k];
}

void SimpleParameterVector::zero()
{
  for (size_t k=0; k<m_buffer.size(); k++) m_buffer[k] = (ParameterScalar)0;
}

ParameterScalar SimpleParameterVector::dot(const ParameterVector *other) const
{
  const SimpleParameterVector* othe=dynamic_cast <const SimpleParameterVector*> (other);
  if (this->variance() * othe->variance() > 0) throw std::logic_error("mismatching co/contravariance");
  ParameterScalar result=0;
  for (size_t k=0; k<m_buffer.size(); k++) result += m_buffer[k] * othe->m_buffer[k];
  return result;
}

void SimpleParameterVector::scal(ParameterScalar a)
{
  for (size_t k=0; k<m_buffer.size(); k++) m_buffer[k] *= a;
}

void SimpleParameterVector::put(ParameterScalar * const buffer, size_t length, size_t offset)
{
    for (size_t k=0; k<length; k++)
        m_buffer[k+offset] = buffer[k];
}

void SimpleParameterVector::get(ParameterScalar *buffer, size_t length, size_t offset) const
{
    for (size_t k=0; k<length; k++)
        buffer[k] = m_buffer[k+offset];
}

size_t SimpleParameterVector::size() const { return this->m_buffer.size();}

std::string SimpleParameterVector::str() const {
    std::ostringstream os; os << "SimpleParameterVector object:";
    for (size_t k=0; k<size(); k++)
        os <<" "<< m_buffer[k];
    os << std::endl;
    return os.str();
}
