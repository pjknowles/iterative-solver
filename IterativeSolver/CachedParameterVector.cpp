#include "CachedParameterVector.h"
#include <stdexcept>
#include <string.h>

using namespace IterativeSolver;

CachedParameterVector::CachedParameterVector(size_t length) :
  m_size(length)
{
    init();
}

CachedParameterVector::CachedParameterVector(const CachedParameterVector& source)
{
    init();
    *this = source;
}

void CachedParameterVector::init()
{
  m_file = new Storage();
  m_cacheMax=m_cacheOffset=m_cacheEmpty=std::numeric_limits<size_t>::max();
  setCacheSize(1024);
}

void CachedParameterVector::setCacheSize(size_t length)
{
  flushCache(true);
  m_cacheSize = length;
  m_cache.resize(m_cacheSize);
  m_cacheOffset=m_cacheEmpty;
  m_cacheDirty=false;
}

CachedParameterVector::~CachedParameterVector() {delete m_file;}


CachedParameterVector& CachedParameterVector::operator=(const CachedParameterVector& other)
{
  std::vector<ParameterScalar> buffer(m_cacheSize);
  m_size=other.m_size;
  if (false) {
  for (size_t block=0; block<m_size; block+=buffer.size()) {
      size_t bs=std::min(buffer.size(),m_size-block);
      other.read(&buffer[0], bs, block);
      write(&buffer[0], bs, block);
    }
    } else
    for (size_t k=0; k<m_size; k++)
      (*this)[k] = other[k];
  setVariance(other.variance());
  return *this;
}

void CachedParameterVector::axpy(ParameterScalar a, const ParameterVector* other)
{
  const CachedParameterVector* othe=dynamic_cast <const CachedParameterVector*> (other);
  if (this->variance() != othe->variance()) throw std::logic_error("mismatching co/contravariance");
  if (false) {
      flushCache();
  std::vector<ParameterScalar> buffer(m_cacheSize);
  std::vector<ParameterScalar> buffero(m_cacheSize);
  for (size_t block=0; block<m_size; block+=buffer.size()) {
      size_t bs=std::min(buffer.size(),m_size-block);
      read(&buffer[0], bs, block);
      othe->read(&buffero[0], bs, block);
      if (true) {
          for (size_t k=0; k<bs; k++) buffer[k] += a*buffero[k];
          write(&buffer[0], bs, block);
        }
      else
        for (size_t k=0; k<bs; k++) (*this)[k+block] += a*buffero[k];
    }
    } else {
        for (size_t k=0; k<m_size; k++) (*this)[k] += a*(*othe)[k];
    }
}

void CachedParameterVector::zero()
{
  if (false) {
      flushCache();
  std::vector<ParameterScalar> buffer(m_cacheSize);
  for (size_t k=0; k<m_cacheSize; k++) buffer[k] += 0;
  for (size_t block=0; block<m_size; block+=buffer.size()) {
      size_t bs=std::min(buffer.size(),m_size-block);
      write(&buffer[0], bs, block);
    }
    } else {
  for (size_t k=0; k<m_size; k++) (*this)[k] = 0;
    }
}


ParameterScalar CachedParameterVector::dot(const ParameterVector *other) const
{
  const CachedParameterVector* othe=dynamic_cast <const CachedParameterVector*> (other);
  if (this->variance() != othe->variance()) throw std::logic_error("mismatching co/contravariance");
  ParameterScalar result=0;
  if (false) {
  std::vector<ParameterScalar> buffer(m_cacheSize);
  std::vector<ParameterScalar> buffero(m_cacheSize);
  for (size_t block=0; block<m_size; block+=buffer.size()) {
      size_t bs=std::min(buffer.size(),m_size-block);
      read(&buffer[0], bs, block);
      othe->read(&buffero[0], bs, block);
      for (size_t k=0; k<bs; k++) result += buffer[k] * buffero[k];
    }
    } else {
      for (size_t k=0; k<m_size; k++)
        result += (*this)[k] * (*othe)[k];
    }
  return result;
}

std::string CachedParameterVector::str() const {
    std::ostringstream os; os << "CachedParameterVector object:";
    flushCache();
//    std::cout << "@ in str, m_cacheDirty="<<m_cacheDirty<<std::endl;
    for (size_t k=0; k<size(); k++) {
//        std::cout << "k="<<k<<", m_cacheDirty="<<m_cacheDirty<<std::endl;
        os <<" "<< (*this)[k];
      }
    os << std::endl;
    return os.str();
}

void CachedParameterVector::put(ParameterScalar * const buffer, size_t length, size_t offset)
{
  flushCache();
  write(buffer,length,offset);
  m_cacheOffset=m_cacheEmpty;
  if (length+offset > m_size) m_size = length+offset;
}

void CachedParameterVector::get(ParameterScalar *buffer, size_t length, size_t offset) const
{
  flushCache();
    read(buffer,length,offset);
}

void CachedParameterVector::flushCache(bool force) const
{
        if (m_cacheOffset==m_cacheEmpty) return;
      if (force || (m_cacheDirty && m_cacheSize < size())) {
//          std::cout << "flush Cache offset="<<m_cacheOffset<<" ,length="<<std::min(m_cacheSize,(size_t)size()-m_cacheOffset)<<std::endl;
//          std::cout << "flush buffer begins "<<m_cache[0]<<std::endl;
          write(&m_cache[0],std::min(m_cacheSize,(size_t)size()-m_cacheOffset),m_cacheOffset);
          m_cacheDirty=false;
        }
}

void CachedParameterVector::write(const ParameterScalar* const buffer, size_t length, size_t address) const
{
  m_file->write((char*) buffer,length*sizeof(ParameterScalar),address*sizeof(ParameterScalar));
}

void CachedParameterVector::read(ParameterScalar* buffer, size_t length, size_t address) const
{
  m_file->read((char*) buffer,length*sizeof(ParameterScalar),address*sizeof(ParameterScalar));
}
