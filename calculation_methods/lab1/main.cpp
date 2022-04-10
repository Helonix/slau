#include <iostream>
#include <chrono>
#include "Matrix.h"
#include "Timer.h"

int main() {
    Timer timer;
    auto begin = std::chrono::steady_clock::now();
    for (int i = 0; i < 50; ++i) {
//        std::cout << i + 1 << " iteration:\n";
        Matrix A = Matrix::get_random_symmetrical_matrix(256);
        Matrix x = Matrix::get_random_vector_column(256);
//    Matrix A = Matrix({{3, -2}, {5, 1}});
//    std::vector<std::vector<double>> vec{{-6}, {3}};
//    Matrix b = Matrix(std::move(vec));
        Matrix b = A * x;
        timer.start();
        A.inverse_by_gauss_jordan_method();
        long long int dur = timer.get_time_in_microseconds();
        std::cout << i << ": " << dur << " microseconds\n";
        timer.start();
        Matrix::solve_by_gauss_with_selection_by_columns(A, b);
        dur = timer.get_time_in_microseconds();
        std::cout << i << ": " << dur << " microseconds\n";
//    std::cout << A.get_condition_number_by_cubic_norm() << " microseconds:" << timer.get_time_in_microseconds() << "\n";
//    std::cout << A << "\n" << "\n" << b << '\n';
//    std::cout << "Attempt\n" << (x == Matrix::solve_by_gauss_with_selection_by_columns(A, b)) << "\n" << (x == Matrix::solve_LUP_by_columns(A.get_LUP_by_columns_decomposition(), b)) << "\n";
//        A.inverse_by_gauss_jordan_method();
////        std::cout << "A:\n" << A << "\n";
////        std::cout << "x:\n" << x << "\n";
////        std::cout << "b:\n" << A * x << "\n";
    }
//    Matrix L = Matrix({{1, 0, 0}, {1.0/3, 1, 0}, {1.0/6, 3.0/4, 1}});
//    Matrix U = Matrix({{-36, 54, -14}, {0, -4, 23.0/3}, {0, 0, 7.0/12}});
//    Matrix P = Matrix({{0, 1, 0}, {0, 0, 1}, {1, 0, 0}});
//    Matrix test = Matrix({{2, 5, 7}, {6, 3, 4}, {5, -2, -3}});
////    Matrix b = Matrix({{-20}, {-4}, {-118}});
//    std::cout << test << "\n" << test.transpose();
    return 0;
}
