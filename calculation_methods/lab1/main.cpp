#include <iostream>
#include "Matrix.h"

int main() {
    for (int i = 0; i < 1; ++i){
        std::cout << i + 1 << " iteration:\n";
        Matrix A = Matrix::get_random_symmetrical_matrix(6);
        Matrix x = Matrix::get_random_vector_column(3);
        A.inverse_by_gauss_jordan_method();
//        std::cout << "A:\n" << A << "\n";
//        std::cout << "x:\n" << x << "\n";
//        std::cout << "b:\n" << A * x << "\n";
    }
    return 0;
}
