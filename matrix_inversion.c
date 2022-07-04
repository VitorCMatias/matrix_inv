#include "matrix_inversion.h"

int LUPDecompose(double **A, int N, double Tol, int *P)
{
    /* INPUT: A - array of pointers to rows of a square matrix having dimension N
     *        Tol - small tolerance number to detect failure when the matrix is near degenerate
     * OUTPUT: Matrix A is changed, it contains a copy of both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
     *        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1
     *        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N,
     *        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S
     */
    int i, j, k, imax;
    double maxA, *ptr, absA;

    for (i = 0; i <= N; i++)
        P[i] = i; // Unit permutation matrix, P[N] initialized with N

    for (i = 0; i < N; i++)
    {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k][i])) > maxA)
            {
                maxA = absA;
                imax = k;
            }

        // printf("AAAAAA %lf\n", maxA);
        if (maxA < Tol)
            return 0; // failure, matrix is degenerate

        if (imax != i)
        {
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

        for (j = i + 1; j < N; j++)
        {
            A[j][i] /= A[i][i];

            for (k = i + 1; k < N; k++)
                A[j][k] -= A[j][i] * A[i][k];
        }
    }

    return 1; // decomposition done
}

void LUPSolve(double **A, int *P, double *b, int N, double *x)
{
    /* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
     * OUTPUT: x - solution vector of A*x=b
     */
    for (int i = 0; i < N; i++)
    {
        x[i] = b[P[i]];

        for (int k = 0; k < i; k++)
            x[i] -= A[i][k] * x[k];
    }

    for (int i = N - 1; i >= 0; i--)
    {
        for (int k = i + 1; k < N; k++)
            x[i] -= A[i][k] * x[k];

        x[i] /= A[i][i];
    }
}

void LUPInvert(double **A, int *P, int N, double **IA)
{
    /* INPUT: A,P filled in LUPDecompose; N - dimension
     * OUTPUT: IA is the inverse of the initial matrix
     */

    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < N; i++)
        {
            IA[i][j] = P[i] == j ? 1.0 : 0.0;

            for (int k = 0; k < i; k++)
                IA[i][j] -= A[i][k] * IA[k][j];
        }

        for (int i = N - 1; i >= 0; i--)
        {
            for (int k = i + 1; k < N; k++)
                IA[i][j] -= A[i][k] * IA[k][j];

            IA[i][j] /= A[i][i];
        }
    }
}

double LUPDeterminant(double **A, int *P, int N)
{
    /* INPUT: A,P filled in LUPDecompose; N - dimension.
     * OUTPUT: Function returns the determinant of the initial matrix
     */

    double det = A[0][0];

    for (int i = 1; i < N; i++)
        det *= A[i][i];

    return (P[N] - N) % 2 == 0 ? det : -det;
}

void print_matrix(double **matix)
{
    for (int i = 0; i < MATRIX_SIZE; ++i)
    {
        for (int j = 0; j < MATRIX_SIZE; ++j)
            printf("%lf\t", matix[i][j]);
        printf("\n");
    }
}

void deallocate_matrix(double **matix)
{
    for (int i = 0; i <= MATRIX_SIZE; ++i)
        free(matix[i]);
    free(matix);
}

void inicializar_matriz_teste(double **matix)
{
    /*
    matix[0][0] = 1.0;
    matix[0][1] = 1.0;
    matix[0][2] = 1.0;
    matix[1][0] = 1.0;
    matix[1][1] = 2.0;
    matix[1][2] = 1.0;
    matix[2][0] = 1.0;
    matix[2][1] = 1.0;
    matix[2][2] = 2.0;*/

    matix[0][0] = -4.0;
    matix[0][1] = -4.0;
    matix[0][2] = -4.0;
    matix[0][3] = 4.0;
    matix[1][0] = -4.0;
    matix[1][1] = -4.0;
    matix[1][2] = 4.0;
    matix[1][3] = -4.0;
    matix[2][0] = -4.0;
    matix[2][1] = 4.0;
    matix[2][2] = -4.0;
    matix[2][3] = -4.0;
    matix[3][0] = 4.0;
    matix[3][1] = -4.0;
    matix[3][2] = -4.0;
    matix[3][3] = -4.0;

}

void allocate_matrix(double **matrix, int size)
{
	for (int i = 0; i <= size; i++)
		matrix[i] = malloc(size * sizeof(double *));
}
