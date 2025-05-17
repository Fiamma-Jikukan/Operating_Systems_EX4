#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>

#define M_PI 3.14159265358979323846
struct timeval start, end;

typedef struct {
    double **result;
    double **m1;
    double **m2;
    int N, L;
    int row, col;
} RowArgs;

double **mult_by_loop(double **m1, double **m2, int M, int N, int L);

double **mult_by_row(double **m1, double **m2, int M, int N, int L);

double **mult_by_element(double **m1, double **m2, int M, int N, int L);

double drand();

double normal();

double **make_matrix(int M, int N);

void print_matrix(double **matrix, int m, int n, char *title);

void set_timer();

void print_elapsed_time();

void free_matrix(double **matrix, int M, int N);

void *row_thread(void *arg);

void *element_thread(void *arg);


int main(int argc, char *argv[]) {
    if (argc < 4) {
        printf("Usage: %s M N L to multiply an MxN matrix by an NxL matrix\n", argv[0]);
        return (EXIT_SUCCESS);
    }
    // create matrices
    int m = atoi(argv[1]);
    int n = atoi(argv[2]);
    int l = atoi(argv[3]);
    double **m1 = make_matrix(m, n);
    print_matrix(m1, m, n, "Matrix 1");
    double **m2 = make_matrix(n, l);
    print_matrix(m2, n, l, "Matrix 2");

    // multiply with a three nested loop
    set_timer();
    double **m3 = mult_by_loop(m1, m2, m, n, l);
    print_elapsed_time();
    print_matrix(m3, m, l, "m1 x m2 by loop");

    // multiply each row with M threads
    set_timer();
    double **m4 = mult_by_row(m1, m2, m, n, l);
    print_elapsed_time();
    print_matrix(m4, m, l, "m1 x m2 by row");

    // multiply each element with M*L threads
    set_timer();
    double **m5 = mult_by_element(m1, m2, m, n, l);
    print_elapsed_time();
    print_matrix(m5, m, l, "m1 x m2 by element");

    free_matrix(m1, m, n);
    free_matrix(m2, n, l);
    free_matrix(m3, m, l);
    free_matrix(m4, m, l);
    free_matrix(m5, m, l);
    return (EXIT_SUCCESS);
}

double **mult_by_loop(double **m1, double **m2, int M, int N, int L) {
    double **result_matrix = (double **) malloc(sizeof(double *) * M);
    for (int i = 0; i < M; i++) {
        result_matrix[i] = malloc(sizeof(double) * L);
    }

    // Compute result[i][j] = sum over k of m1[i][k] * m2[k][j]
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < L; j++) {
            double sum = 0.0;
            for (int k = 0; k < N; k++) {
                sum += m1[i][k] * m2[k][j];
            }
            result_matrix[i][j] = sum;
        }
    }
    return result_matrix;
}

double **mult_by_row(double **m1, double **m2, int M, int N, int L) {
    double **result = malloc(M * sizeof(double *));
    for (int i = 0; i < M; i++)
        result[i] = malloc(L * sizeof(double));

    pthread_t *threads = malloc(M * sizeof(pthread_t)); // make array of threads

    for (int i = 0; i < M; i++) {
        // allocate and fill args for this row
        RowArgs *args = malloc(sizeof(RowArgs));
        args->result = result;
        args->m1 = m1;
        args->m2 = m2;
        args->N = N;
        args->L = L;
        args->row = i;

        if (pthread_create(&threads[i], NULL, row_thread, args) != 0) {
            perror("pthread_create");
            exit(EXIT_FAILURE);
        }
    }

    // wait for all rows to finish
    for (int i = 0; i < M; i++) {
        pthread_join(threads[i], NULL);
    }
    free(threads);
    return result;
}

double **mult_by_element(double **m1, double **m2, int M, int N, int L) {
    double **result = malloc(M * sizeof(double *));
    for (int i = 0; i < M; i++) {
        result[i] = malloc(L * sizeof(double));
    }

    pthread_t *threads = malloc(M * L * sizeof(pthread_t)); // make array of threads
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < L; j++) {
            RowArgs *args = malloc(sizeof(RowArgs));
            args->result = result;
            args->m1 = m1;
            args->m2 = m2;
            args->N = N;
            args->L = L;
            args->row = i;
            args->col = j;
            int idx = i * L + j;
            if (pthread_create(&threads[idx], NULL, element_thread, args) != 0) {
                perror("pthread_create");
                exit(EXIT_FAILURE);
            }
        }
    }
    // wait for all elements to finish
    for (int idx = 0; idx < M * L; idx++) {
        pthread_join(threads[idx], NULL);
    }
    free(threads);
    return result;
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
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i][j] = normal();
        }
    }
    return matrix;
}

void print_matrix(double **matrix, int m, int n, char *title) {
    printf("%s\n", title);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void set_timer() {
    gettimeofday(&start, NULL);
}

void print_elapsed_time() {
    gettimeofday(&end, NULL);
    long time_between_calls = (end.tv_usec - start.tv_usec) + (end.tv_sec - start.tv_sec) * 1000000;
    printf("Time it took to multiply the matrices (in microseconds): %ld\n", time_between_calls);
}

void free_matrix(double **matrix, int M, int N) {
    for (int i = 0; i < M; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

void *row_thread(void *arg) {
    RowArgs *a = (RowArgs *) arg;
    for (int j = 0; j < a->L; j++) {
        double sum = 0.0;
        for (int k = 0; k < a->N; k++) {
            sum += a->m1[a->row][k] * a->m2[k][j];
        }
        a->result[a->row][j] = sum;
    }
    free(a); // clean up arguments struct
    return NULL;
}

void *element_thread(void *arg) {
    RowArgs *a = (RowArgs *) arg;
    double sum = 0.0;
    for (int j = 0; j < a->N; j++) {
        sum += a->m1[a->row][j] * a->m2[j][a->col];
    }
    a->result[a->row][a->col] = sum;
    free(a);
    return NULL;
}
