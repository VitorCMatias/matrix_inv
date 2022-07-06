#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#define MATRIX_SIZE 3

// {{1, 2, 3}, {0, 1, 4}, {0, 0, 1}};
class Matrix {
private:
    int size;
    double **a = nullptr;
    double **Ia = nullptr;
    int *P = nullptr;
    int LU_out = 0;

private:

    int LUPDecompose(double Tol);

    void LUPInvert();

    double LUPDeterminant();

    double **allocate(double **matrix);

    void deallocate(double **matrix);

public:
    Matrix(int size);

    void print(double **A) const;

    void inicializar_matriz_teste();

    void invert();

    void erase();
};
