#ifndef _SOLVER_CPP

#include "Solveur.h"
#include <iostream>

Solver::Solver() {}

Solver::~Solver() {}


void EigenSolver::setSystemMatrix(Eigen::SparseMatrix<double,Eigen::RowMajor> systemMatrix)
{
  Eigen::SparseMatrix<double,Eigen::ColMajor> A(systemMatrix);
  //std::cout<<"Chargement //////////////////"<<std::endl;
  // Compute the ordering permutation vector from the structural pattern of A
  _solver.analyzePattern(A);
  // Compute the numerical factorization
  _solver.factorize(A);
}

Eigen::SparseVector<double> EigenSolver::solve(Eigen::SparseVector<double> RHS)
{
  //Use the factors to solve the linear system
  return _solver.solve(RHS);
}

void BiCGSTABSolver::setSystemMatrix(Eigen::SparseMatrix<double, Eigen::RowMajor> systemMatrix)
{
    _solver.compute(systemMatrix);
}

Eigen::SparseVector<double> BiCGSTABSolver::solve(Eigen::SparseVector<double> RHS)
{
    return _solver.solve(RHS);
}



#define _SOLVER_CPP
#endif