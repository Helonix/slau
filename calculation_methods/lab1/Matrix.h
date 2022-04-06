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
  Colour console_output_colour = NONE;
  unsigned int rows_amount_, columns_amount_;
  std::vector<std::vector<double>> matrix_;


 public:
  Matrix(unsigned int rows_amount, unsigned int columns_amount, double default_value);
  explicit Matrix(unsigned int square_matrix_size, double default_value);
  explicit Matrix(std::vector<std::vector<double>>&& matrix);
  void make_ones_on_main_diag();
  void add(Matrix& lhs);
  void transpose();
  friend std::ostream& operator<< (std::ostream& out, const Matrix& matrix);
  void set_console_text_colour(Colour colour);
  Matrix operator+ (const Matrix& rhs);
};

#endif //CALCULATION_METHODS_LAB1_CMAKE_BUILD_DEBUG_MATRIX_H_
