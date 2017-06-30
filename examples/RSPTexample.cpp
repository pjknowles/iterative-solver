#include "IterativeSolver/ISRSPT.h"
#include "IterativeSolver/PagedParameterVector.h"
#include "IterativeSolver/SimpleParameterVector.h"

using namespace LinearAlgebra;

//  typedef SimpleParameterVector pv;
  typedef PagedParameterVector pv;

// Perturbatively find lowest eigensolution of M(i,j) = alpha*(i+1)*delta(i,j) + i + j
static double n; // dimension of problem
static double alpha; // separation of diagonal elements

static void _matrix_residual(const ParameterVectorSet & psx, ParameterVectorSet & outputs, std::vector<ParameterScalar> shift=std::vector<ParameterScalar>(), bool append=false) {
  std::vector<ParameterScalar> psxk(n);
  std::vector<ParameterScalar> output(n);
  for (size_t k=0; k<psx.size(); k++) {
      psx[k]->get(&(psxk[0]),n,0);
      if (append)
          outputs[k]->get(&output[0],n,0);
      for (size_t i=0; i<n; i++) {
          if (append)
            output[i] += (alpha*(i+1))*psxk[i];
          else
            output[i] = (alpha*(i+1))*psxk[i];
          for (size_t j=0; j<n; j++)
            output[i] += (i+j)*psxk[j];
        }
      outputs[k]->put(&output[0],n,0);
    }
}

static void _matrix_preconditioner(const ParameterVectorSet & psg, ParameterVectorSet & psc, std::vector<ParameterScalar> shift=std::vector<ParameterScalar>(), bool append=false) {
  std::vector<ParameterScalar> psck(n);
  std::vector<ParameterScalar> psgk(n);
  for (size_t k=0; k<psc.size(); k++) {
      psg[k]->get(&psgk[0],n,0);
      if (shift[k]==0)
          for (size_t i=0; i<n; i++)
            psck[i] = psgk[i]*(2*i+alpha*(i+1));
      else if (append) {
          psc[k]->get(&psck[0],n,0);
          for (size_t i=1; i<n; i++)
            psck[i] -= psgk[i]/(shift[k]+2*i+alpha*(i+1));
        } else {
          for (size_t i=1; i<n; i++)
            psck[i] =- psgk[i]/(shift[k]+2*i+alpha*(i+1));
          psck[0] =0;
        }
      psc[k]->put(&psck[0],n,0);
    }
}

int main(int argc, char *argv[])
{
  alpha=100;
  n=10;
  RSPT solver(&_matrix_residual,&_matrix_preconditioner);
  solver.m_verbosity=1;
  solver.m_roots=1;
  ParameterVectorSet g;
  ParameterVectorSet x;
  for (size_t root=0; root<solver.m_roots; root++) {
      x.push_back(std::make_shared<pv>(n));
      g.push_back(std::make_shared<pv>(n));
      x.back()->zero(); double one=1; x.back()->put(&one,1,root); // initial guess
    }
  if (not solver.solve(g,x)) std::cout << "Failure"<<std::endl;
      xout << "Variational eigenvalue "<<solver.eigenvalues().front()<<std::endl;
      for (size_t k=0; k<=solver.iterations(); k++) {
          xout << "E("<<k<<") = "<<solver.incremental_energies()[k]<<", cumulative="<<solver.energy(k)<<", error="<<solver.energy(k)-solver.eigenvalues()[0]<<std::endl;
      }
  std::cout << "Error=";
  for (size_t root=0; root<solver.m_roots; root++)
    std::cout <<solver.errors()[root]<<" ";
  std::cout <<"after "<<solver.iterations()<<" iterations"<<std::endl;
}
