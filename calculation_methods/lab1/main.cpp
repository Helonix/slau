#include <iostream>
#include <fstream>
#include <chrono>
#include "Matrix.h"
#include "Timer.h"

int main() {
    int TEST_MATRICES_AMOUNT = 100;
    int MATRIX_SIZE = 256;
    Timer timer;
    Matrix A, real_x, b, found_x;
    std::vector<std::vector<long long int>> timer_results(TEST_MATRICES_AMOUNT, std::vector<long long int>());
    for (int i = 0; i < TEST_MATRICES_AMOUNT; ++i) {
        A = Matrix::get_random_symmetrical_matrix(MATRIX_SIZE);
        real_x = Matrix::get_random_vector_column(MATRIX_SIZE);
        b = A * real_x;

    }
    std::ofstream out("report.txt");
    return 0;
}
