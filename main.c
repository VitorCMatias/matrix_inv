#include <stdio.h>
#include "matrix_inversion.h"

int main()
{
    double **a = malloc(MATRIX_SIZE * sizeof(double *));
    double **Ia = malloc(MATRIX_SIZE * sizeof(double *));
    int *P = malloc(MATRIX_SIZE + 1 * sizeof(int));
    int LU_out = 0;
    double det = 0;

    for (int i = 0; i < MATRIX_SIZE; ++i)
        a[i] = malloc(MATRIX_SIZE * sizeof(double *));
    for (int i = 0; i < MATRIX_SIZE; ++i)
        Ia[i] = malloc(MATRIX_SIZE * sizeof(double *));

    inicializar_matriz_teste(a);
    print_matrix(a);

    LUPDecompose(a, 3, 1.5, P);
    det = LUPDeterminant(a, P, MATRIX_SIZE);
    LUPInvert(a, P, MATRIX_SIZE, Ia);

    printf("INVERSE\n");
    print_matrix(Ia);

    deallocate_matrix(a);
    deallocate_matrix(Ia);
    free(P);

    return 0;
}
