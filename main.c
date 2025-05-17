#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define M_PI 3.14159265358979323846

double **mult_by_loop(double **m1, double **m2, int M, int N, int L);

double **mult_by_row(double **m1, double **m2, int M, int N, int L);

double **mult_by_element(double **m1, double **m2, int M, int N, int L);

double drand();

double normal();

double **make_matrix(int M, int N);

void print_matrix(double **matrix, int m, int n, char *title);

void set_timer();

void print_elapsed_time();

void free_matrix(double **matrix1, double **matrix2, int M, int N, int L);

int main(int argc, char *argv[]) {
    if (argc < 4) {
        printf("Usage: %s M N L to multiply an MxN matrix by an NxL matrix\n", argv[0]);
        return (EXIT_SUCCESS);
    }
    int m = atoi(argv[1]);
    int n = atoi(argv[2]);
    int l = atoi(argv[3]);
    double **m1 = make_matrix(m, n);
    print_matrix(m1, m, n, "Matrix 1");
    double **m2 = make_matrix(n, l);
    print_matrix(m2, n, l, "Matrix 2");
    set_timer();
    double **m3 = mult_by_loop(m1, m2, m, n, l);
    print_elapsed_time();
    print_matrix(m3, m, l, "m1 x m2 by loop");
    set_timer();
    double **m4 = mult_by_row(m1, m2, m, n, l);
    print_elapsed_time();
    print_matrix(m4, m, l, "m1 x m2 by row");
    set_timer();
    double **m5 = mult_by_element(m1, m2, m, n, l);
    print_elapsed_time();
    print_matrix(m5, m, l, "m1 x m2 by element");
    free_matrix(m1, m2, m, n, l);
    return (EXIT_SUCCESS);
}

double **mult_by_loop(double **m1, double **m2, int M, int N, int L) {
}

double **mult_by_row(double **m1, double **m2, int M, int N, int L) {
}

double **mult_by_element(double **m1, double **m2, int M, int N, int L) {
}

// uniform distribution (0..1)
double drand() {
    return ((random() + 1.0) / (RAND_MAX + 1.0));
}

// normal distribution, centered on 0, std dev 1
double normal() {
    return (sqrt(-2 * log(drand())) * cos(2 * M_PI * drand()));
}

double **make_matrix(int M, int N) {
    double **matrix = (double **) malloc(sizeof(double *) * M);
    for (int i = 0; i < M; i++) {
        matrix[i] = (double *) malloc(sizeof(double) * N);
    }
}

void print_matrix(double **matrix, int m, int n, char *title) {
    printf("%s\n", title);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%p ", &matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void set_timer() {
}

void print_elapsed_time() {
}


void free_matrix(double **matrix1, double **matrix2, int M, int N, int L) {
    for (int i = 0; i < M; i++) {
        free(matrix1[i]);
    }
    free(matrix1);
    for (int i = 0; i < N; i++) {
        free(matrix2[i]);
    }
    free(matrix2);
}
