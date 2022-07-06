#include "matrix_inversion.hpp"

Matrix::Matrix(int size) {
    this->size = size;

    this->a = allocate(this->a);
    this->Ia = allocate(this->Ia);
    this->P = new int[size + 1];
}

int Matrix::LUPDecompose(double Tol) {
    /* INPUT: A - array of pointers to rows of a square matrix having dimension N
     *        Tol - small tolerance number to detect failure when the matrix is near degenerate
     * OUTPUT: Matrix A is changed, it contains a copy of both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
     *        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1
     *        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N,
     *        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S
     */
    int i, j, k, imax;
    double maxA, *ptr, absA;
    int N = this->size;
    int *P = this->P;
    double **A = this->a;

    for (i = 0; i <= N; i++)
        P[i] = i; // Unit permutation matrix, P[N] initialized with N

    for (i = 0; i < N; i++) {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k][i])) > maxA) {
                maxA = absA;
                imax = k;
            }

        // printf("AAAAAA %lf\n", maxA);
        if (maxA < Tol)
            return 0; // failure, matrix is degenerate

        if (imax != i) {
            // pivoting P
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;

            // pivoting rows of A
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;

            // counting pivots starting from N (for determinant)
            P[N]++;
        }

        for (j = i + 1; j < N; j++) {
            A[j][i] /= A[i][i];

            for (k = i + 1; k < N; k++)
                A[j][k] -= A[j][i] * A[i][k];
        }
    }

    return 1; // decomposition done
}

void Matrix::LUPInvert() {
    /* INPUT: A,P filled in LUPDecompose; N - dimension
     * OUTPUT: IA is the inverse of the initial matrix
     */

    double **A = this->a;
    double **IA = this->Ia;
    int *P = this->P;
    int N = this->size;

    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            IA[i][j] = P[i] == j ? 1.0 : 0.0;

            for (int k = 0; k < i; k++)
                IA[i][j] -= A[i][k] * IA[k][j];
        }

        for (int i = N - 1; i >= 0; i--) {
            for (int k = i + 1; k < N; k++)
                IA[i][j] -= A[i][k] * IA[k][j];

            IA[i][j] /= A[i][i];
        }
    }
}

double Matrix::LUPDeterminant() {
    /* INPUT: A,P filled in LUPDecompose; N - dimension.
     * OUTPUT: Function returns the determinant of the initial matrix
     */
    double **A = this->a;
    int *P = this->P;
    int N = this->size;

    double det = A[0][0];

    for (int i = 1; i < N; i++)
        det *= A[i][i];

    return (P[N] - N) % 2 == 0 ? det : -det;
}

void Matrix::print(double **A) const {
    const int size = this->size;

    printf("[");
    for (int i = 0; i < size; ++i) {
        printf("[");
        for (int j = 0; j < size; ++j){
            printf("%lf", A[i][j]);
            j == size - 1 ? printf("]") : printf(",");
        }
        i == size - 1 ? printf("]") : printf(",");
    }
}

void Matrix::inicializar_matriz_teste() {
    double **matrix = this->a;

    matrix[0][0] = 1.0;
    matrix[0][1] = 1.0;
    matrix[0][2] = 1.0;
    matrix[1][0] = 1.0;
    matrix[1][1] = 2.0;
    matrix[1][2] = 1.0;
    matrix[2][0] = 1.0;
    matrix[2][1] = 1.0;
    matrix[2][2] = 2.0;
/*
    matrix[0][0] = -4.0;
    matrix[0][1] = -4.0;
    matrix[0][2] = -4.0;
    matrix[0][3] = 4.0;
    matrix[1][0] = -4.0;
    matrix[1][1] = -4.0;
    matrix[1][2] = 4.0;
    matrix[1][3] = -4.0;
    matrix[2][0] = -4.0;
    matrix[2][1] = 4.0;
    matrix[2][2] = -4.0;
    matrix[2][3] = -4.0;
    matrix[3][0] = 4.0;
    matrix[3][1] = -4.0;
    matrix[3][2] = -4.0;
    matrix[3][3] = -4.0;*/

}

double **Matrix::allocate(double **matrix) {
    const int length = this->size;

    matrix = new double *[length];
    for (int i = 0; i < length; ++i) {
        matrix[i] = new double[length];
    }
    return matrix;
}

void Matrix::invert() {
    double det = 0;

    det = LUPDeterminant();
    if (det != 0) {
        LUPDecompose(1.0);
        LUPInvert();

        printf("INVERSE\n");
        print(this->Ia);
    } else {
        printf("matrix is not inversible, det = 0\r\n");
    }
}

void Matrix::erase() {
    deallocate(a);
    deallocate(Ia);
    delete[]P;
}

void Matrix::deallocate(double **matrix) {

    for (int i = 0; i < this->size; ++i)
        delete matrix[i];
    delete matrix;
}