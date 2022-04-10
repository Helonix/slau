#include <iostream>
#include "Matrix.h"

int main() {
//    for (int i = 0; i < 1; ++i){
//        std::cout << i + 1 << " iteration:\n";
    Matrix A = Matrix::get_random_symmetrical_matrix(10);
    Matrix x = Matrix::get_random_vector_column(10);
//    Matrix A = Matrix({{3, -2}, {5, 1}});
//    std::vector<std::vector<double>> vec{{-6}, {3}};
//    Matrix b = Matrix(std::move(vec));
    Matrix b = A * x;
//    std::cout << A << "\n" << "\n" << b << '\n';
    std::cout << x << "\n" << Matrix::solve_by_gauss_with_selection_by_columns(A, b);
//        A.inverse_by_gauss_jordan_method();
////        std::cout << "A:\n" << A << "\n";
////        std::cout << "x:\n" << x << "\n";
////        std::cout << "b:\n" << A * x << "\n";
//    }
//    Matrix L = Matrix({{1, 0, 0}, {1.0/3, 1, 0}, {1.0/6, 3.0/4, 1}});
//    Matrix U = Matrix({{-36, 54, -14}, {0, -4, 23.0/3}, {0, 0, 7.0/12}});
//    Matrix P = Matrix({{0, 1, 0}, {0, 0, 1}, {1, 0, 0}});
//    Matrix test = Matrix({{2, 5, 7}, {6, 3, 4}, {5, -2, -3}});
////    Matrix b = Matrix({{-20}, {-4}, {-118}});
//    std::cout << test << "\n" << test.transpose();
    return 0;
}
