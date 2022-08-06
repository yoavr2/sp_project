#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "spkmeans.h" /*our library*/

/* check that we free al unnecessary memory */


void free_mat(double **mat, int m){
    int row;
    for (row = 0; row < m; row++){
        free(mat[row]);
    }
    free(mat);
}

double calc_dist(double *x, double *y, int dim){/*calculate (x-y)^2*/


    double tmp, sum;
    int i;
    sum = 0.0;
    tmp = 0.0;
    


    for(i = 0; i<dim; i++){/*calculate (x-centroid)^2 for each coordinate and sum them all*/

        tmp = x[i]-y[i];
        tmp = pow(tmp, 2);
        sum += tmp;

    }/*end of for*/


    return sum;

    /*return -1;*/

}/*end of function calc_dist*/

double calc_exp_euc(double* x, double* y, int N){
    double euc, res;

    euc = calc_dist(x, y, N);
    res = sqrt(fabs(euc));
    res = exp(res*(-1.0/2.0));
    
    return res;
}/*end of function calc_exp_euc*/

double** wam(double** mat, int N){
    
    double** wei_adj_mat;
    int i, j;
    wei_adj_mat = (double **)malloc(N * sizeof(double*));
    if(!wei_adj_mat){/*malloc failed*/
        printf("An Error Has Occurred\n");
        exit(1);
    }

    for(i=0; i<N; i++){
        wei_adj_mat[i] = (double *)malloc(N * sizeof(double));
        if(!wei_adj_mat){/*malloc failed*/
            printf("An Error Has Occurred\n");
            exit(1);
        }
    }

    for(i=0; i<N; i++){
        for(j=0; j<N; j++){
            if (i==j){
                wei_adj_mat[i][j] = 0;
            }
            else if(i > j){
                wei_adj_mat[i][j] = wei_adj_mat[j][i];
            }
            else{
                wei_adj_mat[i][j] = calc_exp_euc(mat[i], mat[j], N);
            }
        }
    }

    return wei_adj_mat;
        
}/*end of function wam*/

double sum_line(double* line, int N){
    int i;
    double sum;
    sum = 0.0;
    
    for(i=0; i<N; i++){
        sum += line[i];
    }

    return sum;
}


double** ddg(double** mat, int N){
    double** dia_deg_mat,** wei_adj_mat;
    int i, j;
    dia_deg_mat = (double **)malloc(N * sizeof(double*));
    if(!dia_deg_mat){/*malloc failed*/
        printf("An Error Has Occurred\n");
        exit(1);
    }

    for(i=0; i<N; i++){
        dia_deg_mat[i] = (double *)malloc(N * sizeof(double));
        if(!dia_deg_mat){/*malloc failed*/
            printf("An Error Has Occurred\n");
            exit(1);
        }
    }

    wei_adj_mat = wam(mat, N);
    
    for(i=0; i<N; i++){
        for(j=0; j<N; j++){
            if (i != j){
                dia_deg_mat[i][j] = 0;
            }
            else{
                dia_deg_mat[i][j] = sum_line(wei_adj_mat[i], N);
            }
        }
    }
    free_mat(wei_adj_mat, N);
    return dia_deg_mat;
}/*end of function ddg*/

double** calc_mat_sqrt(double **dia_deg_mat, int N){
    double **dia_deg_mat_sqrt;
    int i, j;

    dia_deg_mat_sqrt = (double **)malloc(N * sizeof(double*));
    if(!dia_deg_mat_sqrt){/*malloc failed*/
        printf("An Error Has Occurred\n");
        exit(1);
    }

    for(i=0; i<N; i++){
        dia_deg_mat_sqrt[i] = (double *)malloc(N * sizeof(double));
        if(!dia_deg_mat_sqrt){/*malloc failed*/
            printf("An Error Has Occurred\n");
            exit(1);
        }
    }

    for(i=0; i<N; i++){
        for(j=0; j<N; j++){
                if(j!=i){
                    dia_deg_mat_sqrt[i][j] = 0;
                }
                else{/* i == j */
                    dia_deg_mat_sqrt[i][j] = 1.0/sqrt(dia_deg_mat[i][j]);
                }

        }
    }

    return dia_deg_mat_sqrt;
}

double** mat_mult(double** mat1, double** mat2, int N){
    int i,j,k;
    double** rslt;

    rslt = (double **)malloc(N * sizeof(double*));
    if(!rslt){/*malloc failed*/
        printf("An Error Has Occurred\n");
        exit(1);
    }

    for(i=0; i<N; i++){
        rslt[i] = (double *)malloc(N * sizeof(double));
        if(!rslt){/*malloc failed*/
            printf("An Error Has Occurred\n");
            exit(1);
        }
    }


     for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            rslt[i][j] = 0;
 
            for (k = 0; k < N; k++) {
                
                rslt[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
     }
     return rslt;
}

double **mat_mult3(double **mat1, double ** mat2, double **mat3, int N){
    double **res, **mult1;
    
    mult1 = mat_mult(mat1, mat2, N);
    res = mat_mult(mult1, mat3, N);

    return res;

}

double** lnorm(double** data, int N){
    double** dia_deg_mat,** wei_adj_mat,** norm;
    double ** mult1,** mult2,** dia_deg_sqrt_mat;
    int i, j;

    dia_deg_mat = ddg(data, N);
    wei_adj_mat = wam(data, N);
    dia_deg_sqrt_mat = calc_mat_sqrt(dia_deg_mat, N);

    norm = (double **)malloc(N * sizeof(double*));
    if(!norm){/*malloc failed*/
        printf("An Error Has Occurred\n");
        exit(1);
    }

    for(i=0; i<N; i++){
        norm[i] = (double *)malloc(N * sizeof(double));
        if(!norm){/*malloc failed*/
            printf("An Error Has Occurred\n");
            exit(1);
        }
    }

    mult1 = mat_mult(dia_deg_sqrt_mat, wei_adj_mat, N);
    mult2 = mat_mult(mult1, dia_deg_sqrt_mat, N);

    for(i=0; i<N; i++){

            for(j=0; j<N; j++){
                if(i==j){
                    norm[i][j] = 1 - mult2[i][j];
                }
                else{/* j != i */
                    norm[i][j] = 0 - mult2[i][j];
                }
            }
        }

    free(dia_deg_mat);
    free(wei_adj_mat);
    free(mult1);
    free(mult2);
    free(dia_deg_sqrt_mat);

    return norm;
    

}


int sign(double num){
    if (num >= 0){
        return 1;
    }
    return -1;
}

int *find_pivot(double **mat, int N){  /* mat is NxN */
    int i, j, row, column;
    double max_val = mat[0][1];
    int *res = calloc(2, sizeof(int));
    i = 0;
    j = 0;
    for (row = 0; row < N; row++){
        for (column = 0; column < N; column++){
            if (row != column && fabs(mat[row][column]) >= max_val){
                max_val = fabs(mat[row][column]);
                i = row;
                j = column;
            }
        }
    }
    res[0] = i;
    res[1] = j;
    return res;
}


double **Rotation_mat(double **S, int N){  /*returns P the Jacobi rotation matrix of S*/
    double c, s, teta, t, A_ii, A_jj, A_ij;
    double **P;
    int *pivot;
    int row, column, i, j;

    pivot = find_pivot(S, N);

    P = calloc(N, sizeof(double*));
    for (row = 0; row < N; row++){
        P[row] = calloc(N, sizeof(double));
    }

    A_jj = S[pivot[1]][pivot[1]];
    A_ii = S[pivot[0]][pivot[0]];
    A_ij = S[pivot[0]][pivot[1]];
    teta = (A_jj-A_ii)/(2*A_ij);

    t = sign(teta)/(fabs(teta) + sqrt(teta*teta + 1));
    c = 1/sqrt(t*t + 1);
    s = t*c;



    /*i = fmin(pivot[0], pivot[1]);
    j = fmax(pivot[0], pivot[1]);*/
    i = pivot[0];
    j = pivot[1];
    for (row = 0; row < N; row++){
        for (column = 0; column < N; column++){
            if (row == column && row != pivot[0] && row != pivot[1]){ /* case of 1 */
                P[row][column] = 1;
            }
            if (((row == column) && (row == pivot[0]))|| ((row == column) && (row == pivot[1]))){    /* case of c*/
                P[row][column] = c;
            }
            if (row == i && column == N-1-i){  /*case of s  might not be max min*/
                P[row][column] = s;
            }
            if (row == j && column == N-1-j){   /*case of -s    might not be max min*/
                P[row][column] = -1*s;
            }

        }
    }
    return P;

}

double **transpose(double **mat, int N){
    double **transposed;
    int row, column;
    transposed = calloc(N, sizeof(double*));
    for (row = 0; row < N; row++){
        transposed[row] = calloc(N, sizeof(double));
    }
    for (row = 0; row < N; row++){
        for (column = 0; column < N; column++){
            transposed[row][column] = mat[column][row];
        }
    }
    return transposed;

}



double off(double **mat, int N){
    double res;
    int row, column;

    res = 0;

    for (row = 0; row < N; row++){
        for (column = 0; column < N; column++){
            if (row != column){
                res += mat[row][column] * mat[row][column];
            }
        }
    }
    return res;
}

void print_mat(double **mat, int N){
    int row, column;
    for (row = 0; row < N; row++){
        for (column = 0; column < N-1; column++){
            printf("%f, ", mat[row][column]);
        }
        printf("%f\n", mat[row][N-1]);
    }
}



double **Jacoby(double **A, int N){
    double eps, dif, offA, offB;
    double **P, **P_t, **B;
    int counter;


    eps = 1.0 * pow(10, -5);
    counter = 0;
    dif = 10;
    offA = off(A, N);

    while (dif > eps && counter < 100){
        P = Rotation_mat(A, N);
        P_t = transpose(P, N);
        B = mat_mult3(P_t, A, P, N);
        offB = off(B, N);
        dif = fabs(offA - offB);
        free_mat(A, N);
        A = B;
        offA = offB;
        counter += 1;
    }

    
    return A;
    
}

void test_jacoby(){
    double **mat, **A;
    mat = calloc(2, sizeof(double*));
    mat[0] = calloc(2, sizeof(double));
    mat[1] = calloc(2, sizeof(double));
    mat[0][0] = 1;
    mat[0][1] = 2;
    mat[1][0] = 3;
    mat[1][1] = 4;

    A = Jacoby(mat, 2);
    print_mat(A, 2);
    free_mat(A, 2);
    /*free_mat(mat);*/
}

void test_wam(){
    double **mat, **A;
    mat = calloc(2, sizeof(double*));
    mat[0] = calloc(2, sizeof(double));
    mat[1] = calloc(2, sizeof(double));
    mat[0][0] = 1;
    mat[0][1] = 2;
    mat[1][0] = 3;
    mat[1][1] = 4;
    printf("wam:\n");
    A = wam(mat, 2);
    print_mat(A, 2);
    free_mat(A, 2);
    free_mat(mat, 2);
}

void test_ddg(){
    double **mat, **A, **B;
    int i, j;
    mat = calloc(3, sizeof(double*));
    mat[0] = calloc(3, sizeof(double));
    mat[1] = calloc(3, sizeof(double));
    mat[2] = calloc(3, sizeof(double));
    for (i = 0; i < 3; i++){
        for (j = 0; j<3; j++){
            mat[i][j] = i*3+j;
        }
    }
    printf("ddg:\n");
    B = ddg(mat, 3);
    print_mat(B, 3);
    printf("wam:\n");
    A = wam(mat, 3);
    print_mat(A, 3);
    free_mat(B, 3);
}

void test_lnorm(){
    double **mat, **A;
    int i, j;
    mat = calloc(3, sizeof(double*));
    mat[0] = calloc(3, sizeof(double));
    mat[1] = calloc(3, sizeof(double));
    mat[2] = calloc(3, sizeof(double));
    for (i = 0; i < 3; i++){
        for (j = 0; j<3; j++){
            mat[i][j] = i*3+j;
        }
    }

    A = lnorm(mat, 3);
    print_mat(A, 3);
    free_mat(A, 3);
    free_mat(mat, 3);
}


int main(){
    test_lnorm();
    return 0;
}



