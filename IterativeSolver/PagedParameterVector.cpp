#include "PagedParameterVector.h"

#include <stdexcept>
#include <iostream>
#include <string.h>

using namespace IterativeSolver;

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
  char *tmpname=strdup("tmpfileXXXXXX");
  mkstemp(tmpname);
  m_file.open (tmpname, std::ios::out | std::ios::in | std::ios::binary);
  free(tmpname);
}

void PagedParameterVector::setCacheSize(size_t length)
{
    m_cacheSize = length;
    m_cache.resize(m_cacheSize);
    m_cacheOffset=-1;
    m_cacheDirty=false;
}

PagedParameterVector::~PagedParameterVector() {m_file.close();}


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

ParameterScalar PagedParameterVector::at(size_t pos) const
{
    ParameterScalar result;
    read(&result,1,pos);
    return result;
}

ParameterScalar& PagedParameterVector::operator[](size_t pos)
{
  throw std::logic_error("non-const [] cannot be implemented");
}

const ParameterScalar& PagedParameterVector::operator[](size_t pos) const
{
  throw std::logic_error("const [] cannot be implemented");
//  ParameterScalar result;
//  get(&result,1,pos);
//  return result;
}

size_t PagedParameterVector::size() const { return m_size;}

std::string PagedParameterVector::str() const {
    std::ostringstream os; os << "PagedParameterVector object:";
    for (size_t k=0; k<size(); k++)
        os <<" "<< (*this)[k];
    os << std::endl;
    return os.str();
}

void PagedParameterVector::write(ParameterScalar* const buffer, size_t length, size_t address)
{
  m_file.seekg(address*sizeof(ParameterScalar));//,std::ios::beg);
  m_file.write((char*) buffer,length*sizeof(ParameterScalar));
  if (length+address > m_size) m_size = length+address;
}

void PagedParameterVector::put(ParameterScalar * const buffer, size_t length, size_t offset)
{
    write(buffer,length,offset);
}

void PagedParameterVector::read(ParameterScalar* buffer, size_t length, size_t address) const
{
  if (address+length > m_size) throw std::range_error("Storage: attempt to load from beyond end of storage");
  m_file.seekg(address*sizeof(ParameterScalar));//,std::ios::beg);
  m_file.read((char*) buffer,length*sizeof(ParameterScalar));
}

void PagedParameterVector::get(ParameterScalar *buffer, size_t length, size_t offset) const
{
    read(buffer,length,offset);
}
