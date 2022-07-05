#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MATRIX_SIZE 3

// {{1, 2, 3}, {0, 1, 4}, {0, 0, 1}};
class Matrix {
private:
    int size;
    double **a = nullptr;
    double **Ia = nullptr;
    int *P = new int[MATRIX_SIZE + 1];
    int LU_out = 0;

public:
    Matrix(int size);

    int LUPDecompose(double **A, int N, double Tol, int *P);

    void LUPInvert(double **A, int *P, int N, double **IA);

    void LUPSolve(double **A, int *P, double *b, int N, double *x);

    double LUPDeterminant(double **A, int *P, int N);

    void print(double **A);

    void deallocate_matrix(double **matrix);

    void inicializar_matriz_teste();

    double **allocate_matrix(double **matrix);

    void invert();

    void erase();
};
