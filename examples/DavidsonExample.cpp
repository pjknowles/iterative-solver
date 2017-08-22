#include "IterativeSolver/ISDavidson.h"
#include "IterativeSolver/SimpleParameterVector.h"
#include "IterativeSolver/PagedParameterVector.h"
#include "IterativeSolver/CachedParameterVector.h"

using namespace LinearAlgebra;

//  typedef SimpleParameterVector pv;
//  typedef PagedParameterVector pv;
  typedef CachedParameterVector pv;

// Find lowest eigensolutions of M(i,j) = alpha*(i+1)*delta(i,j) + i + j
static size_t n; // dimension of problem
static double alpha; // separation of diagonal elements

struct : IterativeSolverBase::ParameterSetTransformation {
void operator()(const ParameterVectorSet & psx, ParameterVectorSet & outputs, std::vector<scalar> shift=std::vector<scalar>(), bool append=false) const override {
  std::vector<ParameterScalar> psxk(n);
  std::vector<ParameterScalar> output(n);
  for (size_t k=0; k<psx.size(); k++) {
      psx[k]->get(&(psxk[0]),n,0);
      if (append)
          outputs[k]->get(&output[0],n,0);
      else
        for (size_t i=0; i<n; i++) output[i]=0;
      for (size_t i=0; i<n; i++) {
            output[i] += (alpha*(i+1))*psxk[i];
            for (size_t j=0; j<n; j++)
                output[i] += (i+j)*psxk[j];
        }
      outputs[k]->put(&output[0],n,0);
    }
}
} _matrix_residual;

struct : IterativeSolverBase::ParameterSetTransformation {
void operator()(const ParameterVectorSet & psg, ParameterVectorSet & psc, std::vector<scalar> shift=std::vector<scalar>(), bool append=false) const override {
  std::vector<scalar> psck(n);
  std::vector<scalar> psgk(n);
  for (size_t k=0; k<psc.size(); k++) {
      psg[k]->get(&psgk[0],n,0);
      if (append) {
          psc[k]->get(&psck[0],n,0);
          for (size_t i=0; i<n; i++)
            psck[i] -= psgk[i]/(shift[k]+2*i+alpha*(i+1));
        } else {
          for (size_t i=0; i<n; i++)
            psck[i] =- psgk[i]/(shift[k]+2*i+alpha*(i+1));
        }
      psc[k]->put(&psck[0],n,0);
    }
}
} _matrix_preconditioner;

int main(int argc, char *argv[])
{
  alpha=100;
  n=10;
  Davidson solver(_matrix_residual,_matrix_preconditioner);
  solver.m_verbosity=1;
  solver.m_roots=2;
  ParameterVectorSet g;
  ParameterVectorSet x;
  for (size_t root=0; root<solver.m_roots; root++) {
     x.push_back(std::make_shared<pv>(n));
     g.push_back(std::make_shared<pv>(n));
     x.back()->zero(); scalar one=1; x.back()->put(&one,1,root); // initial guess
    }
  if (not solver.solve(g,x)) std::cout << "Failure"<<std::endl;
  std::cout << "Error=";
  for (size_t root=0; root<solver.m_roots; root++)
    std::cout <<solver.errors()[root]<<" ";
  std::cout <<"after "<<solver.iterations()<<" iterations"<<std::endl;
  for (size_t root=0; root<solver.m_roots; root++) {
      std::vector<scalar> buf(n); x[root]->get(&buf[0],n,0);
      std::cout << "Eigenvector:"; for (size_t k=0; k<n; k++) std::cout<<" "<<buf[k]; std::cout<<std::endl;
    }
}
