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
  GAUSS_SOLVE = 1,
  LUP_DECOMPOSITION = 2,
  LUP_SOLVE = 3,
  SQRT_SOLVE = 4,
  RELAXATION_SOLVE = 5
};

#endif //LAB1__STRUCTS_AND_ENUMS_H_
