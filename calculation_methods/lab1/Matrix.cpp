//
// Created by Pavel Emelyanenko on 11.03.2022.
//

#include <iostream>
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
    matrix_ = std::vector(square_matrix_size, std::vector<double>(square_matrix_size, 0));
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


