#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "spkmeans.h" /*our library*/
/* check that we free al unnecessary memory */

/*ENUM*/
enum goal{wam, ddg, lnorm, jacobi};
/*ENUM*/

static int compare (const void * a, const void * b)
{
  if (*(double*)a > *(double*)b) return -1;
  else if (*(double*)a < *(double*)b) return 1;
  else return 0;
}

int find_k(double **mat, int N){
    double *delta, *eigenvalues;
    double max_val;
    int res, i;
    eigenvalues = mat[0];
    qsort(eigenvalues, N, sizeof(double), compare);
    printf("eigenvalues: \n");
    for (i = 0; i < N; i++){
        printf("%f, ", eigenvalues[i]);
    }
    printf("\n");



    delta = calloc(N-1, sizeof(double));
    for (i = 0; i < N-1; i++){
        delta[i] = fabs(eigenvalues[i] - eigenvalues[i+1]);
    }

    max_val = delta[0];
    res = 0;
    for (i = 0; i < N/2; i++){
        if (delta[i] > max_val){
            max_val = delta[i];
            res = i;
        }
    }
    return res + 1;
}

int calcLineLen(char *input_file){
    FILE *fr = fopen(input_file, "r");
    int line_len = 0;
    int max_line = -1;

    if(fr == NULL){
        printf("Invalid Input!\n");
        exit(1);
    }

    while(!(feof(fr))){/*search the maximum line len */

        line_len = 0;

        while(fgetc(fr) != '\n'){

        line_len++;
        if(feof(fr)){
            break;
        }

        }/*end of inner while*/

        if(max_line < line_len){
            max_line = line_len;
        }

        /*if(feof(fr)){
            break;
        }*/

    }/*end of while*/

    max_line++;

    return max_line;

}/*end of function calcLineLen*/


double **new_mat(int N){
    double **mat;
    int i;
    mat = calloc(N, sizeof(double*));
    if(!mat){/*calloc failed*/
        printf("An Error Has Occurred\n");
        exit(1);
    }

    for(i=0; i<N; i++){
        mat[i] = calloc(N, sizeof(double));
        if(!mat){/*calloc failed*/
            printf("An Error Has Occurred\n");
            exit(1);
        }
    }
    return mat;
}

double **I_mat(int N){
    double **mat;
    int i;
    mat = new_mat(N);
    for (i = 0; i < N; i++){
        mat[i][i] = 1;
    }
    return mat;
}

double *get_eigen(double **mat, int N){
    double *eignvalues;
    int i;
    eignvalues = calloc(N, sizeof(double));
    for (i = 0; i < N; i++){
        eignvalues[i] = mat[i][i];
    }
    return eignvalues;
}


double *test_ctopy(){
    double *lst;
    lst = calloc(4, sizeof(double));
    return lst;
}



int *mat_size(char *file_name){
    int *res, rows, columns, i;
    char line[50] = "";
    FILE *fr;
    rows = 0;
    columns = 0;
    fr = fopen(file_name, "r");
    if(fr == NULL){
        printf("Invalid Input!\n");
        exit(1);
    }
    fgets(line, 50, fr);
    rows++;
    i = 0;
    while(line[i] != '\0'){
        if (line[i] == ','){
            columns++;
        }
        i++;
    }
    columns++;


    while (fgets(line, 50, fr) != NULL) {
        rows += 1;
    }

    fclose(fr);

    res = calloc(2, sizeof(int));

    res[0] = rows;
    res[1] = columns;
    printf("rows: %i\n", rows);
    printf("columns: %i\n", columns);
    return res;
}

double **get_mat(char *file_name, int rows, int columns){
    double **mat;
    int i, j;
    int max_chars_in_line;
    char line[50], *token;
    FILE *fr;

    /*max_chars_in_line = calcLineLen( file_name);
    line = (char *)malloc(max_chars_in_line * sizeof(char));*/
    mat = (double **)malloc(rows * sizeof(double*));
    if(!mat){/*malloc failed*/
        printf("An Error Has Occurred\n");
        exit(1);
    }

    for(i=0; i < rows; i++){
        mat[i] = (double *)malloc(columns * sizeof(double));
        if(!mat[i]){/*malloc failed*/
            printf("An Error Has Occurred\n");
            exit(1);
        }
    }

    fr = fopen(file_name, "r");
    if(fr == NULL){
        printf("Invalid Input!\n");
        exit(1);
    }
    i = 0;

    while (fgets(line, 50, fr) != NULL) {
        j = 0;
        token = strtok(line, ",");

        while (token != NULL){
            mat[i][j] = strtod(token, NULL);
            token = strtok(NULL, ",");
            j++;
        }
        i++;

    }

    /*free(line);*/

    return mat;



}

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

double calc_exp_euc(double* x, double* y, int dim){
    double euc, res;

    euc = calc_dist(x, y, dim);
    res = sqrt(fabs(euc));
    res = exp(res*(-1.0/2.0));

    return res;
}/*end of function calc_exp_euc*/

double** wam_func(double **mat, int N, int dim){/*tested on input, worked right*/
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
                wei_adj_mat[i][j] = calc_exp_euc(mat[i], mat[j], dim);
            }
        }
    }

    return wei_adj_mat;

}/*end of function wam_func*/

double sum_line(double* line, int N){
    int i;
    double sum;
    sum = 0.0;

    for(i=0; i<N; i++){
        sum += line[i];
    }

    return sum;
}


double** ddg_func(double **mat, int N, int dim){/*tested on imput, works right*/
    double** dia_deg_mat,** wei_adj_mat;
    int i, j;

    dia_deg_mat = (double **)malloc(N * sizeof(double*));
    if(!dia_deg_mat){
        printf("An Error Has Occurred\n");
        exit(1);
    }

    for(i=0; i<N; i++){
        dia_deg_mat[i] = (double *)malloc(N * sizeof(double));
        if(!dia_deg_mat){
            printf("An Error Has Occurred\n");
            exit(1);
        }
    }

    wei_adj_mat = wam_func(mat, N, dim);

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
}/*end of function ddg_func*/

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

double** lnorm_func(double **mat, int N, int dim){
    double** dia_deg_mat,** wei_adj_mat,** norm;
    double ** mult1,** mult2,** dia_deg_sqrt_mat;
    int i, j;

    dia_deg_mat = ddg_func(mat, N, dim);
    wei_adj_mat = wam_func(mat, N, dim);
    dia_deg_sqrt_mat = calc_mat_sqrt(dia_deg_mat, N);

    norm = (double **)malloc(N * sizeof(double*));
    if(!norm){
        printf("An Error Has Occurred\n");
        exit(1);
    }

    for(i=0; i<N; i++){
        norm[i] = (double *)malloc(N * sizeof(double));
        if(!norm){
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
                else{
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
    j = 1;
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
            if (row == i && column == j){/*N-1-i*/  /*case of s  might not be max min*/
                P[row][column] = s;
            }
            if (row == j && column == i){/*N-1-j*/   /*case of -s    might not be max min*/
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



double **jacobi_func(double **A, int N){/*tested on input, works right*/
    double eps, dif, offA, offB;
    double **P, **P_t, **B, **V, **res;
    double *eigenvalues;
    int counter, i;


    eps = 1.0 * pow(10, -5);
    counter = 0;
    dif = 10;
    offA = off(A, N);
    V = I_mat(N);

    while (dif > eps && counter < 100){
        P = Rotation_mat(A, N);
        V = mat_mult(V, P, N);
        P_t = transpose(P, N);
        B = mat_mult3(P_t, A, P, N);/* it's A' */
        offB = off(B, N);
        dif = fabs(offA - offB);
        free_mat(A, N);
        A = B;
        offA = offB;
        counter += 1;
    }
    eigenvalues = get_eigen(A, N);
    res = calloc(N+1, sizeof(double*));
    for (i = 0; i < N+1; i++){
        res[i] = calloc(N, sizeof(double));
    }
    res[0] = eigenvalues;
    for (i = 0; i < N; i++){
        res[i+1] = V[i];
    }
    free_mat(V, N);

    return res;

}

int heuristic(double **mat, int N, int dim){
    double **lnorm_mat, **jacobi_mat;
    int k;

    lnorm_mat = lnorm_func(mat, N, dim);
    printf("finished lnorm_func in heuristic\n");
    jacobi_mat = jacobi_func(lnorm_mat, N);
    printf("finished jacobi_func in heuristic\n");
    k = find_k(jacobi_mat, N);
    return k;
}

double **spk(double **mat, int N, int dim, int k, int max_iter, int eps){
    double **res;
    printf("k = %i",k);

    res = calloc(2, sizeof(double*));
    res[0] = calloc(2, sizeof(double));
    res[1] = calloc(2, sizeof(double));
    return res;
}


void test_jacoby(){
    double **mat, **A;
    int *size;
    char *file_name;

    file_name = "jacobi_test.txt";
    size = mat_size(file_name);
    mat = get_mat(file_name, size[0], size[1]);

    printf("original mat:\n");
    print_mat(mat, size[0]);

    A = jacobi_func(mat, size[0]);
    printf("Jacobi res:\n");
    print_mat(A, size[0]+1);
}

double **test_pytoc(double **mat, int N){
    print_mat(mat, N);
    //free_mat(mat, N);
    return mat;
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
    printf("wam_func:\n");
    A = wam_func(mat, 2, 2);
    print_mat(A, 2);
    free_mat(A, 2);
    free_mat(mat, 2);
}

/*
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
    printf("ddg_func:\n");
    B = ddg_func(mat, 3);
    print_mat(B, 3);
    printf("wam_func:\n");
    A = wam_func(mat, 3);
    print_mat(A, 3);
    free_mat(B, 3);
}
*/
/*
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

    A = lnorm_func(mat, 3);
    print_mat(A, 3);
    free_mat(A, 3);
    free_mat(mat, 3);
}
*/

/*int numbersOfLines(char *input_file, int line_len){

    int num_of_lines = 0;

    char *buff;

    FILE *fr = fopen(input_file, "r");

    if(fr == NULL){
        printf("Invalid Input!\n");
        exit(1);
    }

    buff = (char *)calloc(line_len+1, sizeof(char));

    if(!buff){
        printf("An Error Has Occurred\n");
        exit(1);
    }

    while(fscanf(fr, "%s", buff ) == 1){/*count the numbers of lines*/

        /*num_of_lines++;

    }/*end of while*/

    /*free(buff);

    return num_of_lines;

    return -1;*/

/*}/*end of function numbersOfLines*/



double ** file_to_mat(char* filename){
    FILE *fr;
    int *row_and_col;
    double **res_mat;
    int i;

    row_and_col = mat_size(filename);

    int num_of_cord = row_and_col[1];
    int data_num = row_and_col[0];

    res_mat = get_mat(filename, data_num, num_of_cord);
    return res_mat;
}



int main(){

    char *file_name;
    int *size, k;
    double **mat, **jacobi_mat, **wam_mat, **ddg_mat, **ddg_sqrt;
    file_name = "tmpFile.txt";
    size = mat_size(file_name);
    mat = get_mat(file_name, size[0], size[1]);

    wam_mat = wam_func(mat, size[0], size[1]);
    ddg_mat = ddg_func(mat, size[0], size[1]);
    ddg_sqrt = calc_mat_sqrt(ddg_mat, size[0]);
    /*k = heuristic(mat, size[0], size[1]);*/

    //printf("k = %d", k);
    /*test_wam();*/

    //jacobi_mat = jacobi_func(mat, size[0]);
    printf("wam_mat: \n");
    print_mat(wam_mat, size[0]);

    printf("\n");

    printf("ddg_mat: \n");
    print_mat(ddg_mat, size[0]);

    printf("\n");

    printf("ddg_sqrt: \n");
    print_mat(ddg_sqrt, size[0]);

    printf("\n");

    test_jacoby();

    return 0;
}



