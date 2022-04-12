#include <iostream>
#include <fstream>
#include "Matrix.h"
#include "Timer.h"


void create_csv_A2plot_file(const std::vector<std::vector<double>>& data){
    {
        std::ofstream out("A2_plot_data.csv");
        out.precision(10);
        out << "0.8,1.0,1.2\n";
        int max_size = std::max((int)data[0].size(), (int)data[1].size());
        max_size = std::max(max_size, (int)data[2].size());
        for (int j = 2; j < max_size; ++j){
            for (int i = 0; i < 3; ++i){
                if (data[i].size() <= j){
                    if (i == 2){
                        out << "\n";
                    }else{
                        out << ",";
                    }
                }else{
                    if (i == 2){
                        out << data[i][j] << "\n";
                    }else{
                        out << data[i][j] << ",";
                    }
                }
            }
        }
    }
}

void create_csv_Aplot_file(const std::vector<std::vector<double>>& data){
    {
        std::ofstream fout("A_plot_data.csv");
        fout.precision(10);
        fout << "0.8,1.0,1.2\n";
        int max_size = std::max((int)data[3].size(), (int)data[4].size());
        max_size = std::max(max_size, (int)data[5].size());
        for (int j = 0; j < max_size; ++j){
            for (int i = 3; i < 6; ++i){
                if (data[i].size() <= j){
                    if (i == 5){
                        fout << "\n";
                    }else{
                        fout << ",";
                    }
                }else{
                    if (i == 5){
                        fout << data[i][j] << "\n";
                    }else{
                        fout << data[i][j] << ",";
                    }
                }
            }
        }
    }
}

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
    Matrix A2;
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
            {1020 - 2 * N, -2924 + 896 * N, 1212 + 9808 * N, -2736 + 98918 * N, 1404 - 11068 * N, -1523 - 8078 * N,
             2625 - 102119 * N, -1327 + 1924 * N};

        Matrix A2_ = Matrix(std::move(vector_A2));
        A2 = A2_.transpose() * A2_;
    }

    Matrix x1 = Matrix::get_random_vector_column(A1.get_rows_amount());
    Matrix x2 = Matrix::get_random_vector_column(A2.get_rows_amount());
    Matrix b1 = A1 * x1;
    Matrix b2 = A2 * x2;

    std::ofstream fout("task9+.txt");
    fout << "Matrix A1:\n" << A1 << "\nMatrix A2:\n" << A2 << "\n";
    long double A1_cond, A2_cond;
    A1_cond = A1.get_condition_number_by_cubic_norm();
    A2_cond = A2.get_condition_number_by_cubic_norm();
    fout << "=============== [Condition Numbers Stats] ===============\n";
    fout << "A1 condition number is " << A1_cond << "\n";
    fout << "A2 condition number is " << A2_cond << "\n\n";

    timer.start();
    A1.inverse_by_gauss_jordan_method();
    duration = timer.get_time_in_microseconds();
    fout << "================ [Inverse Matrix Stats] =================\n";
    fout << "A1 inversion time is " << duration << " mcs\n";
    timer.start();
    A2.inverse_by_gauss_jordan_method();
    duration = timer.get_time_in_microseconds();
    fout << "A2 inversion time is " << duration << " mcs\n\n";

    long double diff;
    timer.start();
    found_x = Matrix::solve_by_gauss_with_selection_by_columns(A1, b1);
    duration = timer.get_time_in_microseconds();
    diff = (x1 - found_x).get_cubic_norm();
    fout << "===================== [Gauss Method] ====================\n";
    fout << "A1 solution difference norm is " << diff << "\n";
    fout << "A1 solution time is " << duration << " mcs\n";

    timer.start();
    found_x = Matrix::solve_by_gauss_with_selection_by_columns(A2, b2);
    duration = timer.get_time_in_microseconds();
    diff = (x2 - found_x).get_cubic_norm();

    fout << "A2 solution difference norm is " << diff << "\n";
    fout << "A2 solution time is " << duration << " mcs\n\n";

    std::pair<Matrix, std::vector<int>> lup;
    timer.start();
    lup = A1.get_LUP_by_columns_decomposition();
    duration = timer.get_time_in_microseconds();

    fout << "============= [LUP-decomposition & solution] ============\n";
    fout << "___-------== ( A1 ) ==-------___\n";
    fout << "LUP-decomposition time is " << duration << " mcs\n";
    fout << "LUP-decomposition:\n";
    Matrix::print_LUP(A1.get_LUP_by_columns_decomposition(), fout);
    timer.start();
    found_x = Matrix::solve_LUP_by_columns(lup, b1);
    duration = timer.get_time_in_microseconds();
    diff = (x1 - found_x).get_cubic_norm();
    fout << "LUP solution time is " << duration << " mcs\n";
    fout << "Solution difference norm is " << diff << "\n\n";

    timer.start();
    lup = A2.get_LUP_by_columns_decomposition();
    duration = timer.get_time_in_microseconds();
    fout << "___-------== ( A2 ) ==-------___\n";
    fout << "LUP-decomposition time is " << duration << " mcs\n";
    fout << "LUP-decomposition:\n";
    Matrix::print_LUP(A2.get_LUP_by_columns_decomposition(), fout);
    timer.start();
    found_x = Matrix::solve_LUP_by_columns(lup, b2);
    duration = timer.get_time_in_microseconds();
    diff = (x2 - found_x).get_cubic_norm();
    fout << "LUP solution time is " << duration << " mcs\n";
    fout << "Solution difference norm is " << diff << "\n\n";

    fout << "================= [LDLT-decomposition] ==================\n";
    fout << "___-------== ( A1 ) ==-------___\n";
    Matrix::print_LDLT(A1.get_LDLT_decomposition(), fout);
    fout << "\n___-------== ( A2 ) ==-------___\n";
    Matrix::print_LDLT(A2.get_LDLT_decomposition(), fout);
    fout << "\n\n";

    bool isSolved = true;
    try {
        timer.start();
        found_x = Matrix::solve_by_sqrt_method(A1, b1);
        duration = timer.get_time_in_microseconds();
        diff = (x1 - found_x).get_cubic_norm();
    } catch (std::runtime_error& err) {
        std::cerr << err.what() << '\n';
        isSolved = false;
        timer.get_time_in_microseconds();
    }
    fout << "================ [Solve by SQRT method] =================\n";
    if (isSolved) {
        fout << "A1 solution difference norm is " << diff << "\n";
        fout << "A1 solution time is " << duration << " mcs\n";
    } else {
        fout << "A1 * x = b1 can't be solved by sqrt method: A1 == A1^T > 0 is not fit\n";
    }

    timer.start();
    found_x = Matrix::solve_by_sqrt_method(A2, b2);
    duration = timer.get_time_in_microseconds();
    diff = (x2 - found_x).get_cubic_norm();

    fout << "A2 solution difference norm is " << diff << "\n";
    fout << "A2 solution time is " << duration << " mcs\n\n";

    std::pair<Matrix, int> relaxation_result;
    timer.start();
    relaxation_result = Matrix::solve_by_relaxation_method(A1, b1, 1 - N / 40, 1e-13);
    duration = timer.get_time_in_microseconds();
    diff = (x1 - relaxation_result.first).get_cubic_norm();

    fout << "============= [Solve by Relaxation method] ==============\n";
    fout << "A1 solution difference norm is " << diff << "\n";
    fout << "A1 solution time is " << duration << " mcs\n";
    fout << "A1 solution iterations number is " << relaxation_result.second << "\n";

    timer.start();
    relaxation_result = Matrix::solve_by_relaxation_method(A2, b2, 1 - N / 40, 1e-6);
    duration = timer.get_time_in_microseconds();
    diff = (x2 - relaxation_result.first).get_cubic_norm();
    fout << "A2 solution difference norm is " << diff << "\n";
    fout << "A2 solution time is " << duration << " mcs\n";
    fout << "A2 solution iterations number is " << relaxation_result.second << "\n\n";


    A = Matrix(matrices[max_cond_index]);
    real_x = Matrix::get_random_vector_column(A.get_rows_amount());
    b = A * real_x;
    Matrix new_b;
    lup = A.get_LUP_by_columns_decomposition();
    long double A_cond = A.get_condition_number_by_cubic_norm();
    long double pr_error, th_error;
    fout << "===== [Vector b Perturbation for Max Cond A Matrix] =====\n";
    double eps = 1e-6;
    for (int i = 1; i < 11; ++i) {
        new_b = b.vector_perturbation(i, eps);
        found_x = Matrix::solve_LUP_by_columns(lup, new_b);
        pr_error = Matrix::get_error(real_x - found_x, real_x);
        th_error = A_cond * Matrix::get_error(b - new_b, b);
        fout << "Iteration " << i << ": \nPerturbation for " << i * eps << "\nPractical error = "
             << pr_error << "\nTheoretical error = " << th_error << "\n";
        fout << "Practical error <= Theoretical error : " << ((th_error - pr_error > 1e-15) ? "True" : "False") << "\n\n";
    }


    A = Matrix(matrices[max_cond_index]);
    real_x = Matrix::get_random_vector_column(A.get_rows_amount());
    b = A * real_x;

    fout << "========= [Relaxation method A2 VS Max Cond A] ==========\n";
    fout << "Diagrams for both matrices are in files A_plot.jpeg and A2_plot.jpeg\n";
    std::vector<std::vector<double>> rel_diff_norms(6, std::vector<double>());


    Matrix::solve_relaxation_and_save_norms(A2, b2, 0.8, 1e-6, rel_diff_norms[0], x2);
    Matrix::solve_relaxation_and_save_norms(A2, b2, 1.0, 1e-6, rel_diff_norms[1], x2);
    Matrix::solve_relaxation_and_save_norms(A2, b2, 1.2, 1e-6, rel_diff_norms[2], x2);
    Matrix::solve_relaxation_and_save_norms(A, b, 0.8, 1e-13, rel_diff_norms[3], real_x);
    Matrix::solve_relaxation_and_save_norms(A, b, 1.0, 1e-13, rel_diff_norms[4], real_x);
    Matrix::solve_relaxation_and_save_norms(A, b, 1.2, 1e-13, rel_diff_norms[5], real_x);

    create_csv_A2plot_file(rel_diff_norms);
    create_csv_Aplot_file(rel_diff_norms);

    system("python ..\\draw_plot.py");


    return 0;
}
