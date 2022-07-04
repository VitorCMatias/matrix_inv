#include <stdio.h>
#include "matrix_inversion.h"

int main()
{
    double **a = malloc(MATRIX_SIZE * sizeof(double *));
    double **Ia = malloc(MATRIX_SIZE * sizeof(double *));
    int *P = malloc(MATRIX_SIZE + 1 * sizeof(int));
    int LU_out = 0;
    double det = 0;

    allocate_matrix(a,MATRIX_SIZE);
    allocate_matrix(Ia,MATRIX_SIZE);

    inicializar_matriz_teste(a);
    print_matrix(a);

    LUPDecompose(a, MATRIX_SIZE, 1.5, P);
    det = LUPDeterminant(a, P, MATRIX_SIZE);
    if(det != 0){
        LUPInvert(a, P, MATRIX_SIZE, Ia);

        printf("INVERSE\n");
        print_matrix(Ia);
    } else{
        printf("matrix is not inversible, det = 0\r\n");
    }


    deallocate_matrix(a);
    deallocate_matrix(Ia);
    free(P);

    return 0;
}
