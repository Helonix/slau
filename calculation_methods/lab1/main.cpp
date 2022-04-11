#include <iostream>
#include <fstream>
#include <chrono>
#include "Matrix.h"
#include "Timer.h"

int main() {
    int TEST_MATRICES_AMOUNT = 100;
    int MATRIX_SIZE = 256;
    double N = 6.0;
    Timer timer;
    Matrix A, real_x, b, found_x;
    std::vector<long long int> timer_results(6);
    std::vector<Matrix> matrices(TEST_MATRICES_AMOUNT);
    std::vector<long double> condition_numbers(TEST_MATRICES_AMOUNT);
    std::vector<int> relaxation_iterations(TEST_MATRICES_AMOUNT);
    std::vector<std::vector<long double>> diff_norms(TEST_MATRICES_AMOUNT, std::vector<long double>(4));
    long long int duration;

    for (int i = 0; i < TEST_MATRICES_AMOUNT; ++i) {
        A = Matrix::get_random_symmetrical_matrix(MATRIX_SIZE);
        real_x = Matrix::get_random_vector_column(MATRIX_SIZE);
        b = A * real_x;
        matrices[i] = A;

        // condition number
        condition_numbers[i] = A.get_condition_number_by_cubic_norm();

        // inverse matrix
        timer.start();
        A.inverse_by_gauss_jordan_method();
        duration = timer.get_time_in_microseconds();
        timer_results[INVERSE_MATRIX] += duration;

        // gauss method solve
        timer.start();
        found_x = Matrix::solve_by_gauss_with_selection_by_columns(A, b);
        duration = timer.get_time_in_microseconds();
        timer_results[GAUSS] += duration;
        diff_norms[i][GAUSS_SOLVE] = (found_x - real_x).get_cubic_norm();

        std::pair<Matrix, std::vector<int>> LUP_decomposition;
        // LUP decomposition
        timer.start();
        LUP_decomposition = A.get_LUP_by_columns_decomposition();
        duration = timer.get_time_in_microseconds();
        timer_results[LUP_DECOMPOSITION] += duration;

        // LUP solve
        timer.start();
        found_x = Matrix::solve_LUP_by_columns(LUP_decomposition, b);
        duration = timer.get_time_in_microseconds();
        timer_results[LUP] += duration;
        diff_norms[i][LUP_SOLVE] = (found_x - real_x).get_cubic_norm();

        // sqrt method solve
        timer.start();
        found_x = Matrix::solve_by_sqrt_method(A, b);
        duration = timer.get_time_in_microseconds();
        timer_results[SQRT] += duration;
        diff_norms[i][SQRT_SOLVE] = (found_x - real_x).get_cubic_norm();

        std::pair<Matrix, int> relaxation_result;
        // relaxation method solve
        timer.start();
        relaxation_result = Matrix::solve_by_relaxation_method(A, b, 1.0 - 6.0 / 40, 1e-15);
        duration = timer.get_time_in_microseconds();
        found_x = relaxation_result.first;
        timer_results[RELAXATION] += duration;
        diff_norms[i][RELAXATION_SOLVE] = (found_x - real_x).get_cubic_norm();
        relaxation_iterations[i] = relaxation_result.second;
    }

    std::vector<long double> average_time(timer_results.size());
    for (int i = 0; i < timer_results.size(); ++i) {
        average_time[i] = static_cast<long double>(timer_results[i]) / TEST_MATRICES_AMOUNT;
    }
    int min_cond_index = 0, max_cond_index = 0;
    long double min_cond = 1e10, max_cond = 0, average_cond = 0;
    for (int i = 0; i < condition_numbers.size(); ++i) {
        average_cond += condition_numbers[i];
        if (condition_numbers[i] - min_cond < 1e-15) {
            min_cond = condition_numbers[i];
            min_cond_index = i;
        }

        if (condition_numbers[i] - max_cond > 1e-15) {
            max_cond = condition_numbers[i];
            max_cond_index = i;
        }
    }
    average_cond /= TEST_MATRICES_AMOUNT;

    int min_iterations = 1e9, max_iterations = 0;
    double average_iterations = 0;

    for (int iterations: relaxation_iterations) {
        average_iterations += iterations;
        if (iterations < min_iterations) {
            min_iterations = iterations;
        }

        if (iterations > max_iterations) {
            max_iterations = iterations;
        }
    }
    average_iterations /= TEST_MATRICES_AMOUNT;

    std::ofstream out("report.txt");
    std::ofstream m_out("min_max_condition_matrices.txt");
    out << "=============== [Condition Numbers Stats] ===============\n";
    out << "Min condition number is " << min_cond << "\n";
    out << "Max condition number is " << max_cond << "\n";
    out << "Average condition number is " << average_cond << "\n\n";
    m_out << "Matrix with minimal condition number:\n" << matrices[min_cond_index]
          << "\nMatrix with maximal condition number:\n" << matrices[max_cond_index];

    out << "================ [Inverse Matrix Stats] =================\n";
    out << "Average time is " << average_time[INVERSE_MATRIX] << " mcs\n\n";

    long double min_diff = 1e10, max_diff = 0, average_diff = 0;
    for (auto diff: diff_norms) {
        average_diff += diff[GAUSS_SOLVE];
        if (diff[GAUSS_SOLVE] < min_diff) {
            min_diff = diff[GAUSS_SOLVE];
        }

        if (diff[GAUSS_SOLVE] > max_diff) {
            max_diff = diff[GAUSS_SOLVE];
        }
    }
    average_diff /= TEST_MATRICES_AMOUNT;

    out << "===================== [Gauss Method] ====================\n";
    out << "Min difference norm is " << min_diff << "\n";
    out << "Max difference norm is " << max_diff << "\n";
    out << "Average difference norm is " << average_diff << "\n";
    out << "Average time is " << average_time[GAUSS] << " mcs\n\n";

    out << "================== [LUP-decomposition] ==================\n";
    out << "Average time is " << average_time[LUP_DECOMPOSITION] << " mcs\n\n";

    min_diff = 1e10;
    max_diff = 0;
    average_diff = 0;
    for (auto diff: diff_norms) {
        average_diff += diff[LUP_SOLVE];
        if (diff[LUP_SOLVE] < min_diff) {
            min_diff = diff[LUP_SOLVE];
        }

        if (diff[LUP_SOLVE] > max_diff) {
            max_diff = diff[LUP_SOLVE];
        }
    }
    average_diff /= TEST_MATRICES_AMOUNT;

    out << "==================== [Solve LUx = b'] ===================\n";
    out << "Min difference norm is " << min_diff << "\n";
    out << "Max difference norm is " << max_diff << "\n";
    out << "Average difference norm is " << average_diff << "\n";
    out << "Average time is " << average_time[LUP] << " mcs\n\n";

    min_diff = 1e10;
    max_diff = 0;
    average_diff = 0;
    for (auto diff: diff_norms) {
        average_diff += diff[SQRT_SOLVE];
        if (diff[SQRT_SOLVE] < min_diff) {
            min_diff = diff[SQRT_SOLVE];
        }

        if (diff[SQRT_SOLVE] > max_diff) {
            max_diff = diff[SQRT_SOLVE];
        }
    }
    average_diff /= TEST_MATRICES_AMOUNT;

    out << "================ [Solve by SQRT method] =================\n";
    out << "Min difference norm is " << min_diff << "\n";
    out << "Max difference norm is " << max_diff << "\n";
    out << "Average difference norm is " << average_diff << "\n";
    out << "Average time is " << average_time[SQRT] << " mcs\n\n";

    min_diff = 1e10;
    max_diff = 0;
    average_diff = 0;
    for (auto diff: diff_norms) {
        average_diff += diff[RELAXATION_SOLVE];
        if (diff[RELAXATION_SOLVE] < min_diff) {
            min_diff = diff[RELAXATION_SOLVE];
        }

        if (diff[RELAXATION_SOLVE] > max_diff) {
            max_diff = diff[RELAXATION_SOLVE];
        }
    }
    average_diff /= TEST_MATRICES_AMOUNT;

    out << "============= [Solve by Relaxation method] ==============\n";
    out << "Min difference norm is " << min_diff << "\n";
    out << "Max difference norm is " << max_diff << "\n";
    out << "Average difference norm is " << average_diff << "\n";
    out << "Average time is " << average_time[RELAXATION] << " mcs\n";
    out << "Min iterations number is " << min_iterations << "\n";
    out << "Max iterations number is " << max_iterations << "\n";
    out << "Average iterations number is " << average_iterations << "\n";

    out.close();
    m_out.close();



    // ============================= [Task 9] ====================================
    Matrix A1;
    {
        std::vector<std::vector<double>> vector_A1 = {{std::pow(N, 2) + 15, N - 1, -1, -2},
                                                      {N - 1, -15 - std::pow(N, 2), -N + 4, -4},
                                                      {-1, -N + 4, std::pow(N, 2) + 8, -N},
                                                      {-2, -4, -N, std::pow(N, 2) + 10}};
        A1 = Matrix(std::move(vector_A1));
    }
    Matrix A2_;
    {
        std::vector<std::vector<double>> vector_A2(8, std::vector<double>(8));
        vector_A2[0][0] = 1;
        for (int i = 1; i < 8; ++i) {
            vector_A2[0][i] = i + N;
        }

        vector_A2[1] = {100 * N, 1000 * N, 10000 * N, 100000 * N, -1000 * N, -10000 * N, -100000 * N, 1};
        vector_A2[2] = {N, -1 + N, -2 + N, -3 + N, -4 + N, -5 + N, -6 + N, -7 + N};
        vector_A2[3] = {N - 1000, 10 * N - 1000, 100 * N - 1000, 1000 * N - 1000, 10000 * N - 1000, -N, -N + 1, -N + 2};
        vector_A2[4] = {-2 * N, 0, -1, -2, -3, -4, -5, -6};
        vector_A2[5] = {N - 2019, -N + 2020, N - 2021, -N + 2022, N - 2023, -N + 2024, N - 2025, -N + 2026};
        vector_A2[6] =
            {2 * N - 2000, 4 * N - 2005, 8 * N - 2010, 16 * N - 2015, 32 * N - 2020, 2019 * N, -2020 * N, 2021 * N};
        vector_A2[7] =
            {1020 - 2 * N, -2924 * 896 * N, 1212 + 9808 * N, -2736 + 98918 * N, 1404 - 11068 * N, -1523 - 8078 * N,
             2625 - 102119 * N, -1327 + 1924 * N};

        A2_ = Matrix(std::move(vector_A2));
    }
    auto A2 = A2_ * A2_.transpose();
    Matrix x1 = Matrix::get_random_vector_column(A1.get_rows_amount());
    Matrix x2 = Matrix::get_random_vector_column(A2.get_rows_amount());
    Matrix b1 = A1 * x1;
    Matrix b2 = A2 * x2;

    std::ofstream fout("task9+.txt");

    return 0;
}
