//
// Created by emely on 10.04.2022.
//

#ifndef LAB1__STRUCTS_AND_ENUMS_H_
#define LAB1__STRUCTS_AND_ENUMS_H_

#include "Matrix.h"
enum Colour {
  GREEN,
  YELLOW,
  RED,
  NONE
};

enum TimeTests {
  INVERSE_MATRIX = 0,
  GAUSS = 1,
  LUP_DECOMPOSITION = 2,
  LUP = 3,
  SQRT = 4,
  RELAXATION = 5
};

enum Norms {
  GAUSS_SOLVE = 0,
  LUP_SOLVE = 1,
  SQRT_SOLVE = 2,
  RELAXATION_SOLVE = 3
};

#endif //LAB1__STRUCTS_AND_ENUMS_H_
