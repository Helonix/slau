#include <iostream>
#include "Matrix.h"

int main() {
    Matrix matrix(5, 0);
    matrix.make_ones_on_main_diag();
    matrix.set_console_text_colour(NONE);
    std::cout << matrix << "\n";
    Matrix matrix1({{1, 2, 3, 4, 5}, {6, 7, 8, 9, 10}, {11, 12, 13, 14, 15}});
    Matrix matrix2({{1, 2, 3, 4, 5}, {6, 7, 8, 9, 10}, {11, 12, 13, 14, 15}});
    std::cout << matrix1 - matrix2;
    return 0;
}
