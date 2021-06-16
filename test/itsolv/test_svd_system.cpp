#include "test.h"
#include <molpro/linalg/itsolv/helper.h>
#include <stdio.h>
#include <molpro/lapacke.h>
#include <vector>
#include <Eigen/Dense>
#include <list>
#include <array>
#include <iostream>

MATCHER_P(FloatNearPointwise, tol, "Out of range") {
    return (std::get<0>(arg)>std::get<1>(arg)-tol && std::get<0>(arg)<std::get<1>(arg)+tol) ;
}

template <typename value_type, int size>
std::vector<value_type> create_test_matrix(){
    std::vector<value_type> m;
    m.resize(size*size, 0);

    for (int i=0; i<size; i++){
        for (int j=0; j<size; j++){
            if( i<=j ){
                float r2 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/1.0));
                m[i+(j*size)] = i == j ? i + 1 : 0.01 * r2 * (i + j);
            }
            else{
                m[i+(j*size)] = m[j+(i*size)];
            }
        }
    }

    return m;
}

void print_matrix( char* desc, lapack_int m, lapack_int n, double* a, lapack_int lda ) {
    lapack_int i, j;
    if (n*lda < 100){
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %8f", a[i*lda+j] );
                printf( "\n" );
        }
    }
}

void print_matrix( char* desc, lapack_int m, lapack_int n, const double* a, lapack_int lda ) {
    lapack_int i, j;
    if (n*lda < 100){
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %8f", a[i*lda+j] );
                printf( "\n" );
        }
    }
}

struct SVDSystem : ::testing::Test {


};


TEST_F(SVDSystem, compare_eigen_solvers){

    // set up test
    const size_t dimension = 5;
    double threshold = 0.5;
    const std::vector<double> m = create_test_matrix<double, dimension>();

    // test lapacke eigensolver
    std::vector<double> eigvecs;
    std::vector<double> eigvals;
    eigvecs.resize(dimension*dimension);
    eigvals.resize(dimension);
    int success = molpro::linalg::itsolv::eigensolver_lapacke_dsyev<double>(m, eigvecs, eigvals, dimension);
    //size_t rank = molpro::linalg::itsolv::get_rank(eigvals, threshold);
    //this causes a linker error for reasons completely unbenknownst to me
    //std::cout << "Rank: " << rank << "\n";

    // reference solution
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> S(m.data(), dimension, dimension);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(S);

    // comparison with SVD
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(S, Eigen::ComputeThinU | Eigen::ComputeThinV);
    svd.setThreshold(threshold);

    // check eigenvalues
    for (int i=0; i<dimension; i++){
        EXPECT_NEAR(eigvals[i], svd.singularValues().reverse()[i], 0.0001 );
        EXPECT_NEAR(es.eigenvalues()[i], eigvals[i], 0.0001 );
        for (int j=0; j<dimension; j++){
            EXPECT_NEAR(abs(eigvecs[i+(j*dimension)]), abs(svd.matrixV().rowwise().reverse()(i,j)), 0.0001 );
            EXPECT_NEAR(abs(eigvecs[i+(j*dimension)]), abs(es.eigenvectors()(i,j)), 0.0001 );
        }
    }

    // print to console if test fails
    if (true){
        char mstr[2] = "m";

        std::cout << "Solution (found using itsolv::eigensolver_lapacke_dsyev) \n";
        char evsstr[100] = "Eigenvalues (computed in lapack)";
        char evcstr[100] = "Eigenvectors";

        print_matrix( mstr, dimension, dimension, m.data(), dimension );
        print_matrix( evsstr, 1, dimension, eigvals.data(), 1 );
        print_matrix( evcstr, dimension, dimension, eigvecs.data(), dimension );

        std::cout << "\n\nSolution (found using Eigen:SelfAdjointEigenSolver): \n";
        std::cout << "Eigenvalues:" << std::endl << es.eigenvalues().transpose() << std::endl << std::endl;
        std::cout << "The matrix of eigenvectors, V:" << std::endl << es.eigenvectors().transpose() << std::endl << std::endl;

        std::cout << "Comparison with SVD: \n";
        std::cout << "singular values (reversed order)\n" << svd.singularValues().transpose().reverse() << std::endl;
        std::cout << "\ndecomposed matrix: (reversed columns) \n" << svd.matrixV().rowwise().reverse().transpose() << std::endl;
        std::cout << "SVD rank " << svd.rank() << " in subspace of dimension " << S.cols() << std::endl;
    }

}