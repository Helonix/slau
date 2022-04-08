//
// Created by emely on 08.04.2022.
//

#ifndef LAB1__SOLAE_H_
#define LAB1__SOLAE_H_

#include "Matrix.h"
class SoLAE {
 private:
  Matrix A_;
  Matrix b_;
  Matrix decomposition_;
  std::vector<int> p_;
 public:
  SoLAE(const Matrix& a, const Matrix& b);
  Matrix solve();
};

#endif //LAB1__SOLAE_H_
