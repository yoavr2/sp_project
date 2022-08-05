#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "spkmeans.h" /*our library*/

/* check that we free al unnecessary memory */


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

double** wam(double** mat){
    
    double** wei_adj_mat;
    int N, i, j;
    N = sizeof(mat)/sizeof(mat[0]);
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


double** ddg(double** mat){
    double** dia_deg_mat,** wei_adj_mat;
    int N, i, j;
    N = sizeof(mat)/sizeof(mat[0]);
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

    wei_adj_mat = wam(mat);
    
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
    free(wei_adj_mat);
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
}

double** lnorm(double** data){
    double** dia_deg_mat,** wei_adj_mat,** norm;
    double ** mult1,** mult2,** dia_deg_sqrt_mat;
    int N, i, j;

    N = sizeof(data)/sizeof(data[0]);
    dia_deg_mat = ddg(data);
    wei_adj_mat = wam(data);
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

int main(int argc, char* argv){
    


}

