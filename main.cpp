#include "matrix_inversion.hpp"

int main()
{

    Matrix A(3);
    for (int i = 0; i < 10000; i++)
    {
        A.set_random(1,10);
        printf("(MATRIX)");
        A.print(A.get_matrix());
        A.invert();
        printf("(INVERSE)");
        A.print(A.get_inverse());
    }

    /*for (int i = 0; i < 100; ++i) {
        A.set_random(1, 99);
        printf("(MATRIX)");
        A.print(A.get_matrix());
        A.invert();
        printf("(INVERSE)");
        A.print(A.get_inverse());
    }*/

    A.erase();

    return 0;
}
