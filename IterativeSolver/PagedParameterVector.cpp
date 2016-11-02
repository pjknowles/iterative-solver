#include "PagedParameterVector.h"
#include <stdexcept>
#include <string.h>

using namespace LinearAlgebra;

PagedParameterVector::PagedParameterVector(size_t length) :
  m_size(length)
{
    init();
}

PagedParameterVector::PagedParameterVector(const PagedParameterVector& source)
{
    init();
    *this = source;
}

void PagedParameterVector::init()
{
  setCacheSize(1024);
  m_file = new Storage();
}

void PagedParameterVector::setCacheSize(size_t length)
{
    m_cacheSize = length;
    m_cache.resize(m_cacheSize);
    m_cacheOffset=-1;
    m_cacheDirty=false;
}

PagedParameterVector::~PagedParameterVector() {delete m_file;}


PagedParameterVector& PagedParameterVector::operator=(const PagedParameterVector& other)
{
  std::vector<ParameterScalar> buffer(m_cacheSize);
  m_size=other.m_size;
  for (size_t block=0; block<m_size; block+=buffer.size()) {
      size_t bs=std::min(buffer.size(),m_size-block);
      other.read(&buffer[0], bs, block);
      write(&buffer[0], bs, block);
    }
  setVariance(other.variance());
  return *this;
}

void PagedParameterVector::axpy(ParameterScalar a, const ParameterVector* other)
{
  const PagedParameterVector* othe=dynamic_cast <const PagedParameterVector*> (other);
  if (this->variance() != othe->variance()) throw std::logic_error("mismatching co/contravariance");
  std::vector<ParameterScalar> buffer(m_cacheSize);
  std::vector<ParameterScalar> buffero(m_cacheSize);
  for (size_t block=0; block<m_size; block+=buffer.size()) {
      size_t bs=std::min(buffer.size(),m_size-block);
      read(&buffer[0], bs, block);
      othe->read(&buffero[0], bs, block);
      for (size_t k=0; k<bs; k++) buffer[k] += a*buffero[k];
      write(&buffer[0], bs, block);
    }
}

void PagedParameterVector::zero()
{
  std::vector<ParameterScalar> buffer(m_cacheSize);
  for (size_t k=0; k<m_cacheSize; k++) buffer[k] += 0;
  for (size_t block=0; block<m_size; block+=buffer.size()) {
      size_t bs=std::min(buffer.size(),m_size-block);
      write(&buffer[0], bs, block);
    }
}


ParameterScalar PagedParameterVector::dot(const ParameterVector *other) const
{
  const PagedParameterVector* othe=dynamic_cast <const PagedParameterVector*> (other);
  if (this->variance() != othe->variance()) throw std::logic_error("mismatching co/contravariance");
  ParameterScalar result=0;
  std::vector<ParameterScalar> buffer(m_cacheSize);
  std::vector<ParameterScalar> buffero(m_cacheSize);
  for (size_t block=0; block<m_size; block+=buffer.size()) {
      size_t bs=std::min(buffer.size(),m_size-block);
      read(&buffer[0], bs, block);
      othe->read(&buffero[0], bs, block);
      for (size_t k=0; k<bs; k++) result += buffer[k] * buffero[k];
    }
  return result;
}

void PagedParameterVector::scal(ParameterScalar a)
{
  std::vector<ParameterScalar> buffer(m_cacheSize);
  for (size_t block=0; block<m_size; block+=buffer.size()) {
      size_t bs=std::min(buffer.size(),m_size-block);
      read(&buffer[0], bs, block);
      for (size_t k=0; k<bs; k++) buffer[k] *= a;
      write(&buffer[0], bs, block);
    }
}

size_t PagedParameterVector::size() const { return m_size;}

std::string PagedParameterVector::str() const {
    std::ostringstream os; os << "PagedParameterVector object:";
    std::vector<ParameterScalar> buffer(m_cacheSize);
    for (size_t block=0; block<m_size; block+=buffer.size()) {
        size_t bs=std::min(buffer.size(),m_size-block);
        read(&buffer[0], bs, block);
        for (size_t k=0; k<size(); k++)
          os <<" "<< buffer[k];
      }
    os << std::endl;
    return os.str();
}

void PagedParameterVector::put(ParameterScalar * const buffer, size_t length, size_t offset)
{
    write(buffer,length,offset);
}

void PagedParameterVector::get(ParameterScalar *buffer, size_t length, size_t offset) const
{
    read(buffer,length,offset);
}

void PagedParameterVector::write(ParameterScalar* const buffer, size_t length, size_t address)
{
  m_file->write((char*) buffer,length*sizeof(ParameterScalar),address*sizeof(ParameterScalar));
  if (length+address > m_size) m_size = length+address;
}

void PagedParameterVector::read(ParameterScalar* buffer, size_t length, size_t address) const
{
  m_file->read((char*) buffer,length*sizeof(ParameterScalar),address*sizeof(ParameterScalar));
}
