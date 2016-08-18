#include "diiscxx.h"
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>

diiscxx::diiscxx(size_t AmpLength, size_t ResLength, size_t OtherLength, double threshold,DiisMode_type DiisMode)
  : AmpLength_(AmpLength), ResLength_(ResLength), OtherLength_(OtherLength), threshold_(threshold), DiisMode_(DiisMode)
{
  if (ResLength_==0) ResLength_=AmpLength_;
  TotalLength_=AmpLength+ResLength+OtherLength;
  setDiisMode();
}

void diiscxx::setDiisMode(DiisMode_type DiisMode) { DiisMode_ = DiisMode; }

// default I/O implementation
inline void diiscxx::dumpInit()
{
  char *tmpname=strdup("tmpfileXXXXXX");
  mkstemp(tmpname);
  dumpFile.open (tmpname, std::ios::out | std::ios::in | std::ios::binary);
  free(tmpname);
}
inline void diiscxx::dumpEnd()
{
  dumpFile.close();
}

void diiscxx::dump(const double *buffer, size_t length, size_t address)
{
  dumpFile.seekg(address,std::ios::beg);
  dumpFile.write((char*) buffer,length);
}

void diiscxx::load(double *buffer, int index, size_t length, size_t address)
{
  dumpFile.seekg(address,std::ios::beg);
  dumpFile.read((char*) buffer,length);
}

