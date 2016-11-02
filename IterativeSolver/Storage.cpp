#include "Storage.h"
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <stdexcept>

using namespace LinearAlgebra;

#include <unistd.h>
Storage::Storage(size_t lengthHint, int option)
{
  char *tmpname=strdup("tmpfileXXXXXX");
  mkstemp(tmpname);
  m_file.open (tmpname, std::ios::out | std::ios::in | std::ios::binary);
  unlink(tmpname);
  free(tmpname);
  size_=0;
}

Storage::~Storage()
{
  m_file.close();
}

void Storage::write(const char *buffer, size_t length, size_t address)
{
  m_file.seekg(address);//,std::ios::beg);
  m_file.write((char*) buffer,length);
  if (length+address > size_) size_ = length+address;
}

void Storage::read(char *buffer, size_t length, size_t address) const
{
  if (address+length > size_) throw std::range_error("Storage: attempt to load from beyond end of storage");
  m_file.seekg(address);//,std::ios::beg);
  m_file.read((char*) buffer,length);
}

size_t Storage::size() const
{
  return size_;
}
