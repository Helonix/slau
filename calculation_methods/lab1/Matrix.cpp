//
// Created by Pavel Emelyanenko on 11.03.2022.
//

#include <iostream>
#include <iomanip>
#include <random>
#include "Matrix.h"

std::random_device rd;
std::default_random_engine engine(rd());
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
            out << std::setw(matrix.output_precision_ * 2) << std::fixed << matrix.matrix_[i][j];
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

//Matrix::Matrix(const Matrix& other) {
//    this->matrix_ = std::vector(other.matrix_);
//    this->rows_amount_ = other.rows_amount_;
//    this->columns_amount_ = other.columns_amount_;
//}

// TODO: change implementation: no unary minus usage
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

    for (int i = 0; i < A.rows_amount_ - 1; ++i) {
        if (std::abs(A.matrix_[i][i]) < eps_) {
            std::cerr << "Can't find inverse matrix, determinant is equal to 0\n";
            return {this->rows_amount_, this->columns_amount_, 0};
        }
        for (int k = i + 1; k < A.rows_amount_; ++k) {
            if (std::abs(A.matrix_[k][i]) < eps_) {
                continue;
            }
            double ratio = -(A.matrix_[k][i] / A.matrix_[i][i]);
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
        if (std::abs(A.matrix_[i][i]) < eps_) {
            std::cerr << "Can't find inverse matrix, determinant is equal to 0\n";
            return {this->rows_amount_, this->columns_amount_, 0};
        }
        for (int k = i - 1; k >= 0; --k) {
            if (std::abs(A.matrix_[k][i]) < eps_) {
                continue;
            }
            double ratio = -(A.matrix_[k][i] / A.matrix_[i][i]);
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
        A.matrix_[i][i] = 1;
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

// не работает, проблема в индексах, скорее всего
Matrix Matrix::solve_by_gauss_with_selection_by_columns(Matrix A, Matrix b) {
    Matrix x = Matrix(A.rows_amount_, 1, 0);
    std::vector<int> p(A.rows_amount_);

    for (int i = 0; i < A.rows_amount_; ++i) {
        p[i] = i;
    }
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
            std::swap(p[i], p[max_elem_index]);
        }
        for (int k = i + 1; k < A.rows_amount_; ++k) {
            if (std::abs(A.matrix_[k][i]) < A.eps_) {
                continue;
            }
            double ratio = -(A.matrix_[k][i] / A.matrix_[i][i]);
//            A.matrix_[k][i] = -ratio;

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

    return std::pair<Matrix, std::vector<int>>();
}

Matrix::Matrix() {
    this->matrix_ = std::vector<std::vector<double>>(1, std::vector<double>(1, 0));
    this->rows_amount_ = 1;
    this->columns_amount_ = 1;
}
