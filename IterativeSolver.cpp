#include "IterativeSolver.h"
#include <algorithm>

using namespace IterativeSolver;

IterativeSolverBase::IterativeSolverBase(ParameterSetTransformation updateFunction, ParameterSetTransformation residualFunction)
:  updateFunction_(updateFunction), residualFunction_(residualFunction), verbosity_(0), thresh_(1e-6)
{}

IterativeSolverBase::~IterativeSolverBase()
{
}


#include <limits>
//void IterativeSolverBase::addPreconditioner(double *d, double shift, bool absolute)
//{
//    double shifter = shift;
//    if (not absolute) {
//        shifter=std::numeric_limits<double>::min();
//        for (size_t k=0; k<length_; k++)
//          if (d[k]<-shifter) shifter=-d[k];
//        shifter +=  shift;
//        if (shift == 0) shift+=1e-14; // to avoid division by zero
//      }
//    preconditioner_store_ = new Storage(sizeof(double)*length_);
//  std::vector<double> buffer(std::min(buffer_size_,length_));
//  for (size_t block=0; block<length_; block+=buffer.size()) {
//      for (size_t k=0; k<std::min(buffer_size_,length_-block); k++)
//        buffer[k] = 1/(d[block+k] + shift);
//      preconditioner_store_->write(&buffer[0],std::min(buffer.size(),length_-block)*sizeof(double),(block)*sizeof(double));
//        }
//}

//void IterativeSolverBase::update(const ParameterVectorSet &residual, ParameterVectorSet &solution)
//{
//    updateFunction_(residual,solution,0);
//}

bool IterativeSolverBase::iterate(const ParameterVectorSet &residual, ParameterVectorSet &solution)
{
    residuals_.push_back(residual);
    solutions_.push_back(solution);
    updateFunction_(residual,solution,0);
    return calculateError(residual,solution) < thresh_;
}

std::vector<double> IterativeSolverBase::calculateErrors(const ParameterVectorSet &residual, const ParameterVectorSet &solution)
{
    std::vector<double> result;
    for (size_t k=0; k<residual.size(); k++)
        result.push_back(residual.active[k] ? residual[k]*solution[k] : 0);
    return result;
}

double IterativeSolverBase::calculateError(const ParameterVectorSet &residual, const ParameterVectorSet &solution)
{
    std::vector<double> errs=calculateErrors(residual,solution);
    return std::sqrt(std::inner_product(errs.begin(),errs.end(),errs.begin(),0));
}



// default I/O implementation
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
Storage::Storage(size_t lengthHint, int option)
{
  char *tmpname=strdup("tmpfileXXXXXX");
  mkstemp(tmpname);
  dumpFile_.open (tmpname, std::ios::out | std::ios::in | std::ios::binary);
  free(tmpname);
  size_=0;
}

Storage::~Storage()
{
  dumpFile_.close();
}

void Storage::write(const double *buffer, size_t length, size_t address)
{
  dumpFile_.seekg(address);//,std::ios::beg);
  dumpFile_.write((char*) buffer,length);
//  std::cout << "Storage::write at address "<< address << ", length="<<length<<","<< buffer[0]<<std::endl;
  if (length+address > size_) size_ = length+address;
}

void Storage::read(double *buffer, size_t length, size_t address)
{
  if (address+length > size_) throw std::range_error("Storage: attempt to load from beyond end of storage");
  dumpFile_.seekg(address);//,std::ios::beg);
  dumpFile_.read((char*) buffer,length);
//  std::cout << "Storage::read at address "<< address << ", length="<<length<<","<< buffer[0]<<std::endl;
}

size_t Storage::size()
{
  return size_;
}
