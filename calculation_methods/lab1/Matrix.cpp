//
// Created by Pavel Emelyanenko on 11.03.2022.
//

#include <iostream>
#include <iomanip>
#include "Matrix.h"
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

    for (int i = 0; i < rows_amount_; ++i){
        for (int j = 0; j < columns_amount_; ++j){
            matrix_[i][j] += lhs.matrix_[i][j];
        }
    }
}

// TODO: transposition is not finished
void Matrix::transpose() {
    unsigned int new_rows = columns_amount_;
    unsigned int new_columns = rows_amount_;
    if (rows_amount_ != columns_amount_){
        unsigned int new_size = std::max(columns_amount_, rows_amount_);
        matrix_.resize(new_size);
        for (int i = 0; i < matrix_.size(); ++i){
            matrix_[i].resize(new_size);
        }
    }
}

std::ostream& operator<<(std::ostream& out, const Matrix& matrix) {
//    out.precision(7);
    std::ios state(nullptr);
    state.copyfmt(out);
    out.precision(8);
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
            out  << std::setw(20) << std::fixed <<  matrix.matrix_[i][j];
        }
        out << "\n";
    }
    out << "\x1b[0m";
    out.copyfmt(state);
    return out;
}
void Matrix::set_console_text_colour(Colour colour) {
    console_output_colour = colour;
}
Matrix Matrix::operator+(const Matrix& rhs) {
    if (this->rows_amount_ != rhs.rows_amount_ || this->columns_amount_ == rhs.columns_amount_){
        std::cerr << "Inconsistent matrix for operation +\n";
        return *this;
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
    matrix_ = std::vector(matrix);
    rows_amount_ = matrix_.size();
    if (rows_amount_ > 0){
        columns_amount_ = matrix_[0].size();
    }else{
        std::cerr << "Empty matrix\n";
        return;
    }
}
