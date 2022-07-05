#include <stdio.h>
#include "matrix_inversion.hpp"


int main() {
    double **a = nullptr;
    double **Ia = nullptr;
    int *P = new int[MATRIX_SIZE + 1];
    int LU_out = 0;
    double det = 0;

    a = allocate_matrix(a);
    Ia = allocate_matrix(Ia);

    inicializar_matriz_teste(a);
    print_matrix(a);

    LUPDecompose(a, MATRIX_SIZE, 1.5, P);
    det = LUPDeterminant(a, P, MATRIX_SIZE);
    if (det != 0) {
        LUPInvert(a, P, MATRIX_SIZE, Ia);

        printf("INVERSE\n");
        print_matrix(Ia);
    } else {
        printf("matrix is not inversible, det = 0\r\n");
    }

    deallocate_matrix(a);
    deallocate_matrix(Ia);
    delete[]P;

    return 0;
}


