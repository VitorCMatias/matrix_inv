#include "matrix_inversion.hpp"


int main() {

    Matrix A(3);

    A.inicializar_matriz_teste();
    A.invert();
    A.erase();



    return 0;
}


