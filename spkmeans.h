
#ifndef SPKMEANS_H
#define SPKMEANS_H

double *test_ctopy();
double** wam_func(double **mat, int N, int dim);
double** ddg_func(double **mat, int N, int dim);
double** lnorm_func(double **mat, int N, int dim);
double **jacobi_func(double **A, int N);
double **spk(double **mat, int N, int dim, int k, int max_iter, int eps);
int heuristic(double **mat, int N, int dim);
#endif
