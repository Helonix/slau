//
// Created by Pavel Emelyanenko on 11.03.2022.
//

#ifndef LAB1__MATRIX_H_
#define LAB1__MATRIX_H_

#include <vector>
#include <iostream>
#include <random>
#include "structs_and_enums.h"

class Matrix {
 private:
  double eps_ = 1e-7;
  Colour console_output_colour = NONE;
  unsigned int rows_amount_, columns_amount_;
  std::vector<std::vector<double>> matrix_;
  int output_precision_ = 7;
//  std::pair<Matrix, Matrix> last_decomposition_;

 public:
  // Constructors
  Matrix(unsigned int rows_amount, unsigned int columns_amount, double default_value);
  explicit Matrix(unsigned int square_matrix_size, double default_value);
  explicit Matrix(std::vector<std::vector<double>>&& matrix);
  Matrix();
//  Matrix(const Matrix& other);

  // Operator overloading
  Matrix operator+(const Matrix& rhs) const;
  Matrix operator-(const Matrix& rhs) const;
  Matrix operator*(const Matrix& rhs) const;
  Matrix operator-() const;
  bool operator==(const Matrix& rhs) const;
  bool operator!=(const Matrix& rhs) const;
  friend std::ostream& operator<<(std::ostream& out, const Matrix& matrix);

  // Getters
  long double get_cubic_norm() const;
  unsigned int get_rows_amount() const;
  unsigned int get_columns_amount() const;
  long double get_condition_number_by_cubic_norm();
  std::pair<Matrix, std::vector<int>> get_LUP_by_columns_decomposition();

  // Setters
  void set_console_text_colour(Colour colour);
  void set_output_precision(int precision);

  // 1st + 2nd tasks: create matrix & vector
  static Matrix get_random_symmetrical_matrix(unsigned int size);
  static Matrix get_random_vector_column(unsigned int size);

  // 3rd task: inverse matrix (cubic norm and condition number are in getters)
  Matrix inverse_by_gauss_jordan_method();

  // 4th task: gauss with selection by columns
  static Matrix solve_by_gauss_with_selection_by_columns(Matrix A, Matrix b);

  // 5th task: solve LUP by columns
  static Matrix solve_LUP_by_columns(const std::pair<Matrix, std::vector<int>>& lup_decomposition, const Matrix& b);

  // Other methods
//  static std::string LUP_by_columns_to_string(std::pair<Matrix, std::vector<int>> lup);
  static std::string LDLT_to_string(const Matrix& decomposition);
  void make_ones_on_main_diag();
  void add(Matrix& lhs);
  Matrix transpose();
};

#endif //LAB1__MATRIX_H_
