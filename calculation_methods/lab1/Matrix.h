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
  double eps_ = std::numeric_limits<double>::epsilon() * 10;
  Colour console_output_colour = NONE;
  unsigned int rows_amount_, columns_amount_;
  std::vector<std::vector<double>> matrix_;
  int output_precision_ = 12;
  bool save_diff_norms = true;

 public:
  // Constructors
  Matrix(unsigned int rows_amount, unsigned int columns_amount, double default_value);
  explicit Matrix(unsigned int square_matrix_size, double default_value);
  explicit Matrix(std::vector<std::vector<double>>&& matrix);
  Matrix();
  Matrix(const Matrix& other);

  // Operator overloading
  Matrix operator+(const Matrix& rhs) const;
  Matrix operator-(const Matrix& rhs) const;
  Matrix operator*(const Matrix& rhs) const;
  Matrix operator*(int rhs) const;
  Matrix operator-() const;
  bool operator==(const Matrix& rhs) const;
  bool operator!=(const Matrix& rhs) const;
  friend std::ostream& operator<<(std::ostream& out, const Matrix& matrix);

  // Getters
  double get_cubic_norm() const;
  unsigned int get_rows_amount() const;
  unsigned int get_columns_amount() const;
  long double get_condition_number_by_cubic_norm();
  std::pair<Matrix, std::vector<int>> get_LUP_by_columns_decomposition();
  Matrix get_GGT_decomposition();
  Matrix get_LDLT_decomposition();
  static Matrix get_identity_matrix(unsigned int size);

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

  // 6th task: solve by square root method
  static Matrix solve_by_sqrt_method(Matrix A, Matrix b);

  // 7th task: relaxation method
  static std::pair<Matrix, int> solve_by_relaxation_method(Matrix A, Matrix b, double w, double eps);
  static std::pair<Matrix, int> solve_relaxation_and_save_norms(Matrix A, Matrix b, double w, double eps, std::vector<double>& norms_data, const Matrix& real_x);

  // Other methods
  static void print_LUP(const std::pair<Matrix, std::vector<int>>& lup, std::ostream& out);
  static void print_LDLT(const Matrix& decomposition, std::ostream& out);
  void make_ones_on_main_diag();
  void add(Matrix& lhs);
  Matrix transpose();
  Matrix vector_perturbation(int iteration, double value);
  static double get_error(const Matrix& difference, const Matrix& base);
  bool is_positive();
};

#endif //LAB1__MATRIX_H_
