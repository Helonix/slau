//
// Created by emely on 08.04.2022.
//

#include "SoLAE.h"
SoLAE::SoLAE(const Matrix& a, const Matrix& b) : A_(a), b_(b), decomposition_(a) {
    p_ = std::vector<int>(a.columns_amount_);
    for (int i = 0; i < a.columns_amount_; ++i) {
        p_[i] = i + 1;
    }
}
Matrix SoLAE::solve() {
    return Matrix(0, 0, 0);
}
