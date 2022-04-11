//
// Created by Pavel Emelyanenko on 11.03.2022.
//

#include <iostream>
#include <iomanip>
#include "Matrix.h"

std::default_random_engine engine(time(nullptr));
std::uniform_real_distribution<double> dist(-std::pow(2.0, 1.5), std::pow(2.0, 1.5));

/// Common matrix constructor
/// \param rows_amount
/// \param columns_amount
/// \param default_value - some default double value that is used to fill matrix
Matrix::Matrix(unsigned int rows_amount, unsigned int columns_amount, double default_value) {
    matrix_ = std::vector(rows_amount, std::vector<double>(columns_amount, default_value));
    rows_amount_ = rows_amount;
    columns_amount_ = columns_amount;
}

/// Constructor for square matrix
/// \param square_matrix_size
/// \param default_value - some default double value that is used to fill matrix
Matrix::Matrix(unsigned int square_matrix_size, double default_value) {
    matrix_ = std::vector(square_matrix_size, std::vector<double>(square_matrix_size, default_value));
    rows_amount_ = square_matrix_size;
    columns_amount_ = square_matrix_size;
}


void Matrix::make_ones_on_main_diag() {
    if (this->columns_amount_ != this->rows_amount_){
        std::cerr << "Matrix is not square, use another method" << "\n";
        return;
    }

    for (int i = 0; i < rows_amount_; ++i){
        matrix_[i][i] = 1;
    }
}

/// Matrix addition, changes are made for the current matrix
/// \param lhs - matrix that we need to add to the current matrix
void Matrix::add(Matrix& lhs) {
    if (columns_amount_ != lhs.columns_amount_ || rows_amount_ != lhs.rows_amount_){
        std::cerr << "Sizes of matrices are not equal\n";
        return;
    }

    for (int i = 0; i < rows_amount_; ++i) {
        for (int j = 0; j < columns_amount_; ++j) {
            matrix_[i][j] += lhs.matrix_[i][j];
        }
    }
}

Matrix Matrix::transpose() {
    Matrix transposed = Matrix(this->columns_amount_, this->rows_amount_, 0);

    for (int i = 0; i < this->rows_amount_; ++i) {
        for (int j = 0; j < this->columns_amount_; ++j) {
            transposed.matrix_[j][i] = this->matrix_[i][j];
        }
    }

    return transposed;
}

std::ostream& operator<<(std::ostream& out, const Matrix& matrix) {
    std::ios state(nullptr);
    state.copyfmt(out);
    out.precision(matrix.output_precision_);
    switch (matrix.console_output_colour) {
        case GREEN:
            out << "\x1b[32;1;3m";
            break;
        case RED:
            out << "\x1b[31;1;3m";
            break;
        case YELLOW:
            out << "\x1b[33;1;3m";
            break;
        case NONE:
            break;
    }
    for (int i = 0; i < matrix.rows_amount_; ++i){
        for (int j = 0; j < matrix.columns_amount_; ++j){
            out << std::setw(matrix.output_precision_ * 2 + 5) << std::fixed << matrix.matrix_[i][j];
        }
        out << "\n";
    }
    if (matrix.console_output_colour != NONE){
        out << "\x1b[0m";
    }
    out.copyfmt(state);
    return out;
}

void Matrix::set_console_text_colour(Colour colour) {
    console_output_colour = colour;
}

Matrix Matrix::operator+(const Matrix& rhs) const {
    if (this->rows_amount_ != rhs.rows_amount_ || this->columns_amount_ != rhs.columns_amount_){
        std::cerr << "Inconsistent matrix for operation +\n";
        return {this->rows_amount_, this->columns_amount_, 0};
    }
    Matrix result = Matrix(this->rows_amount_, this->columns_amount_, 0);
    for (int i = 0; i < rows_amount_; ++i){
        for (int j = 0; j < columns_amount_; ++j){
            result.matrix_[i][j] = this->matrix_[i][j] + rhs.matrix_[i][j];
        }
    }
    return result;
}

Matrix::Matrix(std::vector<std::vector<double>>&& matrix) {
    matrix_ = std::vector<std::vector<double>>(matrix);
    rows_amount_ = matrix_.size();
    if (rows_amount_ > 0){
        columns_amount_ = matrix_[0].size();
    }else{
        std::cerr << "Empty matrix\n";
        return;
    }
}

Matrix Matrix::operator*(const Matrix& rhs) const {
    if (this->columns_amount_ != rhs.rows_amount_){
        std::cerr << "Matrices are inconsistent\n";
        return {this->rows_amount_, this->columns_amount_, 0};
    }

    Matrix result(this->rows_amount_, rhs.columns_amount_, 0);
    for (int i = 0; i < this->rows_amount_; ++i){
        for (int j = 0; j < rhs.columns_amount_; ++j){
            for (int k = 0; k < this->columns_amount_; ++k){
                result.matrix_[i][j] += this->matrix_[i][k] * rhs.matrix_[k][j];
            }
        }
    }

    return result;
}


Matrix Matrix::operator-() const {
    Matrix result(this->rows_amount_, this->columns_amount_, 0);

    for (int i = 0; i < this->rows_amount_; ++i) {
        for (int j = 0; j < this->columns_amount_; ++j) {
            result.matrix_[i][j] = -this->matrix_[i][j];
        }
    }

    return result;
}

Matrix::Matrix(const Matrix& other) {
    this->matrix_ = std::vector(other.matrix_);
    this->rows_amount_ = other.rows_amount_;
    this->columns_amount_ = other.columns_amount_;
}

Matrix Matrix::operator-(const Matrix& rhs) const {
    return *this + (-rhs);
}

Matrix Matrix::get_random_symmetrical_matrix(unsigned int size) {
    Matrix result(size, size, 0);
    result.make_ones_on_main_diag();

    for (int i = 0; i < size; ++i) {
        for (int j = i + 1; j < size; ++j) {
            result.matrix_[i][j] = dist(engine);
            result.matrix_[j][i] = result.matrix_[i][j];

            result.matrix_[i][i] += std::abs(result.matrix_[i][j]);
            result.matrix_[j][j] += std::abs(result.matrix_[i][j]);
        }
    }

    return result;
}
Matrix Matrix::get_random_vector_column(unsigned int size) {
    Matrix vector(size, 1, 0);

    for (int i = 0; i < size; ++i){
        vector.matrix_[i][0] = dist(engine);
    }

    return vector;
}

long double Matrix::get_cubic_norm() const {
    long double max_row_sum = 0, current_sum = 0;
    for (int i = 0; i < this->rows_amount_; ++i) {
        for (int j = 0; j < this->columns_amount_; ++j) {
            current_sum += std::abs(this->matrix_[i][j]);
        }
        if (current_sum - max_row_sum > eps_) {
            max_row_sum = current_sum;
        }
        current_sum = 0;
    }
    return max_row_sum;
}

// TODO: протестировать дополнительно
Matrix Matrix::inverse_by_gauss_jordan_method() {
    if (this->rows_amount_ != this->columns_amount_) {
        std::cerr << "Matrix must be square\n";
        return {this->rows_amount_, this->columns_amount_, 0};
    }

    Matrix inverse = {this->rows_amount_, this->columns_amount_, 0};
    inverse.make_ones_on_main_diag();
    Matrix A = Matrix(*this);

    // forward Gauss
    double ratio;
    for (int i = 0; i < A.rows_amount_ - 1; ++i) {
        if (std::abs(A.matrix_[i][i]) < eps_) {
            std::cerr << "Can't find inverse matrix, determinant is equal to 0\n";
            return {this->rows_amount_, this->columns_amount_, 0};
        }
        for (int k = i + 1; k < A.rows_amount_; ++k) {
            if (std::abs(A.matrix_[k][i]) < eps_) {
                continue;
            }
            ratio = -(A.matrix_[k][i] / A.matrix_[i][i]);
            for (int j = 0; j < A.columns_amount_; ++j) {
                inverse.matrix_[k][j] += inverse.matrix_[i][j] * ratio;
                if (j >= i) {
                    A.matrix_[k][j] += A.matrix_[i][j] * ratio;
                }
            }
        }
    }

    // Reverse Gauss


    for (int i = A.rows_amount_ - 1; i > 0; --i) {
        for (int k = i - 1; k >= 0; --k) {
            if (std::abs(A.matrix_[k][i]) < eps_) {
                continue;
            }
            ratio = -(A.matrix_[k][i] / A.matrix_[i][i]);
            A.matrix_[k][i] += A.matrix_[i][i] * ratio;
            for (int j = 0; j < A.columns_amount_; ++j) {
                inverse.matrix_[k][j] += inverse.matrix_[i][j] * ratio;
            }
        }
    }

    for (int i = 0; i < A.rows_amount_; ++i) {
        for (int j = 0; j < A.columns_amount_; ++j) {
            inverse.matrix_[i][j] /= A.matrix_[i][i];
        }
    }

    return inverse;
}
unsigned int Matrix::get_rows_amount() const {
    return rows_amount_;
}

unsigned int Matrix::get_columns_amount() const {
    return columns_amount_;
}

long double Matrix::get_condition_number_by_cubic_norm() {
    return this->inverse_by_gauss_jordan_method().get_cubic_norm() * this->get_cubic_norm();
}

void Matrix::set_output_precision(int precision) {
    output_precision_ = precision;
}

Matrix Matrix::solve_by_gauss_with_selection_by_columns(Matrix A, Matrix b) {
    Matrix x = Matrix(A.rows_amount_, 1, 0);

    int max_elem_index = 0;
    double max_elem = 0;
    for (int i = 0; i < A.rows_amount_; ++i) {
        for (int m = i; m < A.rows_amount_; ++m) {
            if (std::abs(A.matrix_[m][i]) - max_elem_index > A.eps_) {
                max_elem = std::abs(A.matrix_[m][i]);
                max_elem_index = m;
            }
        }
        if (max_elem_index != i) {
            std::swap(A.matrix_[i], A.matrix_[max_elem_index]);
            std::swap(b.matrix_[i], b.matrix_[max_elem_index]);
        }
        for (int k = i + 1; k < A.rows_amount_; ++k) {
            if (std::abs(A.matrix_[k][i]) < A.eps_) {
                continue;
            }
            double ratio = -(A.matrix_[k][i] / A.matrix_[i][i]);

            for (int j = i + 1; j < A.columns_amount_; ++j) {
                A.matrix_[k][j] += ratio * A.matrix_[i][j];
            }
            b.matrix_[k][0] += b.matrix_[i][0] * ratio;
        }
    }

    double summary = 0;
    for (int i = A.rows_amount_ - 1; i >= 0; --i) {
        summary = b.matrix_[i][0];
        for (int j = i + 1; j < A.rows_amount_; ++j) {
            summary -= A.matrix_[i][j] * x.matrix_[j][0];
        }
        x.matrix_[i][0] = summary / A.matrix_[i][i];
    }

    return x;
}

std::pair<Matrix, std::vector<int>> Matrix::get_LUP_by_columns_decomposition() {
    std::vector<int> p(this->rows_amount_);

    for (int i = 0; i < this->rows_amount_; ++i) {
        p[i] = i;
    }

    Matrix A = Matrix(*this);

    int max_elem_index = 0;
    double max_elem = 0;
    for (int i = 0; i < A.rows_amount_; ++i) {
        for (int m = i; m < A.rows_amount_; ++m) {
            if (std::abs(A.matrix_[m][i]) - max_elem_index > A.eps_) {
                max_elem = std::abs(A.matrix_[m][i]);
                max_elem_index = m;
            }
        }
        if (max_elem_index != i) {
            std::swap(A.matrix_[i], A.matrix_[max_elem_index]);
            std::swap(p[i], p[max_elem_index]);
        }
        for (int k = i + 1; k < A.rows_amount_; ++k) {
            if (std::abs(A.matrix_[k][i]) < A.eps_) {
                continue;
            }
            double ratio = -(A.matrix_[k][i] / A.matrix_[i][i]);
            A.matrix_[k][i] = -ratio;

            for (int j = i + 1; j < A.columns_amount_; ++j) {
                A.matrix_[k][j] += ratio * A.matrix_[i][j];
            }
        }
    }

    return std::pair<Matrix, std::vector<int>>{A, p};
}

Matrix::Matrix() {
    this->matrix_ = std::vector<std::vector<double>>(1, std::vector<double>(1, 0));
    this->rows_amount_ = 1;
    this->columns_amount_ = 1;
}

bool Matrix::operator==(const Matrix& rhs) const {
    if (columns_amount_ != rhs.columns_amount_ || rows_amount_ != rhs.rows_amount_) {
        return false;
    }

    for (int i = 0; i < this->rows_amount_; ++i) {
        for (int j = 0; j < this->columns_amount_; ++j) {
            if (std::abs(this->matrix_[i][j] - rhs.matrix_[i][j]) > this->eps_) {
                return false;
            }
        }
    }

    return true;
}

bool Matrix::operator!=(const Matrix& rhs) const {
    return !(rhs == *this);
}

Matrix Matrix::solve_LUP_by_columns(const std::pair<Matrix, std::vector<int>>& lup_decomposition, const Matrix& b) {
    Matrix x = Matrix(lup_decomposition.first.rows_amount_, 1, 0);
    Matrix p_b = Matrix(b.rows_amount_, 1, 0);
    Matrix A = Matrix(lup_decomposition.first);
    for (int i = 0; i < p_b.rows_amount_; ++i) {
        p_b.matrix_[i][0] = b.matrix_[lup_decomposition.second[i]][0];
    }

    Matrix y = Matrix(x.rows_amount_, x.columns_amount_, 0);
    double summary;
    for (int i = 0; i < A.rows_amount_; ++i) {
        summary = p_b.matrix_[i][0];
        for (int j = 0; j < i; ++j) {
            summary -= A.matrix_[i][j] * y.matrix_[j][0];
        }
        y.matrix_[i][0] = summary;
    }

    for (int i = A.rows_amount_ - 1; i >= 0; --i) {
        summary = y.matrix_[i][0];
        for (int j = i + 1; j < A.columns_amount_; ++j) {
            summary -= A.matrix_[i][j] * x.matrix_[j][0];
        }
        x.matrix_[i][0] = summary / A.matrix_[i][i];
    }

    return x;
}


std::pair<Matrix, int> Matrix::solve_by_relaxation_method(Matrix A, Matrix b, double w, double eps = 1e-13) {
    Matrix previous_x = Matrix(A.rows_amount_, 1, 1);
    Matrix x = Matrix(A.rows_amount_, 1, 0);
    double ratio;
    int MAX_ITERATIONS = 500;
    int iterations = 0;
    while (iterations < MAX_ITERATIONS && (x - previous_x).get_cubic_norm() > eps) {
        previous_x = Matrix(x);
        for (int i = 0; i < A.rows_amount_; ++i) {
            ratio = w / A.matrix_[i][i];
            x.matrix_[i][0] = (1 - w) * previous_x.matrix_[i][0] + ratio * b.matrix_[i][0];
            for (int j = 0; j < x.rows_amount_; ++j) {
                if (j > i) {
                    x.matrix_[i][0] -= ratio * A.matrix_[i][j] * previous_x.matrix_[j][0];
                } else if (j < i) {
                    x.matrix_[i][0] -= ratio * A.matrix_[i][j] * x.matrix_[j][0];
                }
            }
        }
        ++iterations;
    }

    return std::pair<Matrix, int>{x, iterations};
}

Matrix Matrix::operator*(int rhs) const {
    Matrix result = Matrix(*this);
    for (int i = 0; i < result.rows_amount_; ++i) {
        for (int j = 0; j < result.columns_amount_; ++j) {
            result.matrix_[i][j] *= rhs;
        }
    }
    return result;
}
Matrix Matrix::get_GGT_decomposition() {
    Matrix A(*this);
    double ratio;
    for (int i = 0; i < A.rows_amount_; ++i) {
        for (int k = i + 1; k < A.rows_amount_; ++k) {
            if (std::abs(A.matrix_[k][i]) < A.eps_) {
                continue;
            }
            ratio = -(A.matrix_[k][i] / A.matrix_[i][i]);
            for (int j = i + 1; j < A.columns_amount_; ++j) {
                A.matrix_[k][j] += ratio * A.matrix_[i][j];
            }
        }
    }

    for (int i = 0; i < A.rows_amount_; ++i) {
        ratio = std::sqrt(std::abs(A.matrix_[i][i]));
        A.matrix_[i][i] /= ratio;
        for (int j = i + 1; j < A.columns_amount_; ++j) {
            A.matrix_[i][j] /= ratio;
            A.matrix_[j][i] = A.matrix_[i][j];
        }
    }

    return A;
}

Matrix Matrix::get_LDLT_decomposition() {
    Matrix A(*this);
    double ratio;
    for (int i = 0; i < A.rows_amount_; ++i) {
        for (int k = i + 1; k < A.rows_amount_; ++k) {
            if (std::abs(A.matrix_[k][i]) < A.eps_) {
                continue;
            }
            ratio = -(A.matrix_[k][i] / A.matrix_[i][i]);
            A.matrix_[k][i] = -ratio;

            for (int j = i + 1; j < A.columns_amount_; ++j) {
                A.matrix_[k][j] += ratio * A.matrix_[i][j];
            }
        }
    }

    for (int i = 0; i < A.rows_amount_; ++i) {
        for (int j = i + 1; j < A.columns_amount_; ++j) {
            A.matrix_[i][j] = A.matrix_[j][i];
        }
    }

    return A;
}

Matrix Matrix::solve_by_sqrt_method(Matrix A, Matrix b) {

    Matrix x = Matrix(b.rows_amount_, 1, 0);
    Matrix y = Matrix(b.rows_amount_, 1, 0);

    A = A.get_GGT_decomposition();

    for (int i = 0; i < A.rows_amount_; ++i) {
        y.matrix_[i][0] = b.matrix_[i][0];
        for (int j = 0; j < i; ++j) {
            y.matrix_[i][0] -= A.matrix_[i][j] * y.matrix_[j][0];
        }
        y.matrix_[i][0] /= A.matrix_[i][i];
    }

    for (int i = A.rows_amount_ - 1; i >= 0; --i) {
        x.matrix_[i][0] = y.matrix_[i][0];
        for (int j = i + 1; j < A.columns_amount_; ++j) {
            x.matrix_[i][0] -= A.matrix_[i][j] * x.matrix_[j][0];
        }
        x.matrix_[i][0] /= A.matrix_[i][i];
    }

    return x;
}
void Matrix::print_LDLT(const Matrix& decomposition, std::ostream& out) {
    Matrix L = Matrix(decomposition.rows_amount_, 0);
    Matrix D = Matrix(decomposition.rows_amount_, 0);
    for (int i = 0; i < decomposition.rows_amount_; ++i) {
        D.matrix_[i][i] = decomposition.matrix_[i][i];
        L.matrix_[i][i] = 1;
        for (int j = 0; j < i; ++j) {
            L.matrix_[i][j] = decomposition.matrix_[i][j];
        }
    }

    out << "L matrix:\n" << L << "\nD matrix:\n" << D << "\nL^T matrix:\n" << L.transpose() << "\n";
    out << L * D * L.transpose();
}


