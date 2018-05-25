#include "IterativeSolver.h"
#include "PagedVector.h"
// Find lowest eigensolutions of M(i,j) = alpha*(i+1)*delta(i,j) + i + j
// Storage of vectors distributed and out of memory via PagedVector class
using scalar = double;
using pv = LinearAlgebra::PagedVector<scalar>;
constexpr size_t n=100; // dimension of problem
constexpr double alpha=100; // separation of diagonal elements

void action(const LinearAlgebra::vectorSet<scalar> & psx, LinearAlgebra::vectorSet<scalar> & outputs) {
 std::vector<scalar> psxk(n);
 std::vector<scalar> output(n);
 for (size_t k=0; k<psx.size(); k++) {
  psx[k]->get(&(psxk[0]),n,0);
  for (size_t i=0; i<n; i++) {
   output[i] = (alpha*(i+1))*psxk[i];
   for (size_t j=0; j<n; j++)
    output[i] += (i+j)*psxk[j];
  }
  outputs[k]->put(&output[0],n,0);
 }
}

void update(LinearAlgebra::vectorSet<scalar> & psc, const LinearAlgebra::vectorSet<scalar> & psg, std::vector<scalar> shift=std::vector<scalar>()) {
 std::vector<scalar> psck(n);
 std::vector<scalar> psgk(n);
 for (size_t k=0; k<psc.size(); k++) {
  psg[k]->get(&psgk[0],n,0);
  for (size_t i=0; i<n; i++)
   psck[i] =- psgk[i]/(shift[k]+2*i+alpha*(i+1));
  psc[k]->put(&psck[0],n,0);
 }
}

int main(int argc, char *argv[])
{
  LinearAlgebra::LinearEigensystem<scalar> solver;
  solver.m_verbosity=1;
  solver.m_roots=4;
  solver.m_thresh=1e-6;
  LinearAlgebra::vectorSet<double> g;
  LinearAlgebra::vectorSet<double> x;
  for (int root=0; root<solver.m_roots; root++) {
     x.push_back(std::make_shared<pv>(n));
     g.push_back(std::make_shared<pv>(n));
     x.back()->zero(); scalar one=1; x.back()->put(&one,1,root); // initial guess
    }
  for (auto iter=0; iter<1000; iter++) {
   action(x,g);
   solver.addVector(x,g);
   update(x,g,solver.eigenvalues());
   if (solver.endIteration(x,g)) break;
  }
  std::cout << "Error={ ";
  for (int root=0; root<solver.m_roots; root++)
    std::cout <<solver.errors()[root]<<" ";
  std::cout <<"} after "<<solver.iterations()<<" iterations"<<std::endl;
  for (int root=0; root<solver.m_roots; root++) {
      std::vector<scalar> buf(n); x[root]->get(&buf[0],n,0);
      std::cout << "Eigenvector:"; for (size_t k=0; k<n; k++) std::cout<<" "<<buf[k]; std::cout<<std::endl;
    }
}
