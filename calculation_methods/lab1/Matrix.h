//
// Created by Pavel Emelyanenko on 11.03.2022.
//

#ifndef CALCULATION_METHODS_LAB1_CMAKE_BUILD_DEBUG_MATRIX_H_
#define CALCULATION_METHODS_LAB1_CMAKE_BUILD_DEBUG_MATRIX_H_

#include <vector>

enum Colour{
  GREEN,
  YELLOW,
  RED,
  NONE
};

class Matrix {
 private:
  double eps_ = 1e-7;
  Colour console_output_colour = NONE;
  unsigned int rows_amount_, columns_amount_;
  std::vector<std::vector<double>> matrix_;


 public:
  // Constructors
  Matrix(unsigned int rows_amount, unsigned int columns_amount, double default_value);
  explicit Matrix(unsigned int square_matrix_size, double default_value);
  explicit Matrix(std::vector<std::vector<double>>&& matrix);
  Matrix(const Matrix& other);

  // Operator overloading
  Matrix operator+ (const Matrix& rhs) const;
  Matrix operator- (const Matrix& rhs) const;
  Matrix operator* (const Matrix& rhs) const;
  Matrix operator- () const;
  friend std::ostream& operator<< (std::ostream& out, const Matrix& matrix);

  // Getters
  long double get_cubic_norm() const;

  // Setters
  void set_console_text_colour(Colour colour);

  // 1st task: create matrix & vector
  static Matrix get_random_symmetrical_matrix(unsigned int size);
  static Matrix get_random_vector_column(unsigned int size);

  // Other methods
  void make_ones_on_main_diag();
  void add(Matrix& lhs);
  void transpose();
  Matrix inverse_by_gauss_jordan_method();
};

#endif //CALCULATION_METHODS_LAB1_CMAKE_BUILD_DEBUG_MATRIX_H_
