#pragma once
#include <armadillo>

arma::vec analytic_beam_eigenvalues(double diagonal, double offDiagonal, int n);
arma::vec armadillo_eigenvalues(int n, double dia, double offDia, arma::vec addDia);
int* offDiag(arma::mat A, int* p, int* q, int n);
arma::mat Rotate(arma::mat A, arma::mat R, int k, int l, int n);
arma::mat setupTriDiag(int n, double diagonal, double offDiagonal, arma::vec addDiagonal);
arma::vec JacobiSolver(double tol, int n, int maxiter, double dia, double offDia, arma::vec addDia);
